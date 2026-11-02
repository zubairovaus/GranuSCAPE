"""Cell refinement, nucleus associations and explicit classification support."""

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from .features import measure_profile
from .expansion import compute
from .classification import fit_model, predict
from .marginal import marginal_membership

FEATURES = ["log_soma_area", "log_nucleus_area", "log1p_AF488", "log_DAPI"]


def measure_cells(labels, nuclei, raw, pixel_size):
    boxes, nuclear_boxes = ndi.find_objects(labels), ndi.find_objects(nuclei)
    pad = int(np.ceil(12 / pixel_size))
    records, ranks = [], []
    border = set(
        np.unique(np.r_[labels[0], labels[-1], labels[:, 0], labels[:, -1]])
    ) - {0}
    for sid, box in enumerate(boxes, 1):
        if box is None:
            continue
        y0, y1 = max(0, box[0].start - pad), min(labels.shape[0], box[0].stop + pad)
        x0, x1 = max(0, box[1].start - pad), min(labels.shape[1], box[1].stop + pad)
        for nid in np.unique(nuclei[box][labels[box] == sid]):
            if nid:
                nb = nuclear_boxes[int(nid) - 1]
                y0, y1 = min(y0, nb[0].start), max(y1, nb[0].stop)
                x0, x1 = min(x0, nb[1].start), max(x1, nb[1].stop)
        sl = np.s_[y0:y1, x0:x1]
        profile, rank = measure_profile(
            labels[sl], nuclei[sl], raw["af488"][sl], raw["dapi"][sl], pixel_size, sid
        )
        profile.update(soma_id=sid, touches_frame=sid in border)
        profile["soma_centroid_x_px"] += x0
        profile["soma_centroid_y_px"] += y0
        profile["nuclear_fraction"] = (
            1 - profile["cytoplasm_area_um2"] / profile["soma_area_um2"]
        )
        # Conservative diagnostic flag; not a validated biological classifier.
        profile["mask_review_required"] = profile["nuclear_fraction"] >= 0.8
        profile["measurement_valid"] = not (
            sid in border or profile["mask_review_required"]
        )
        records.append(profile)
        ranks.append(rank)
    return pd.DataFrame(records), ranks


def refine_cells(raw, masks, pixel_size):
    """Evaluate every existing soma seed, independently of area, contrast or MCL membership."""
    selected = np.unique(masks["soma"])
    selected = selected[selected > 0].astype(int).tolist()
    refined, report = compute(
        raw["af488"], masks["soma"], masks["nuclei"], pixel_size, selected
    )
    for record in report:
        sid = record["soma_id"]
        if record.get("review_flags"):
            refined[(refined == sid) & (masks["soma"] != sid)] = 0
            record["applied"] = False
        else:
            record["applied"] = record["status"] in (
                "expanded_requires_review",
                "unchanged",
            )
    if not np.array_equal(refined[masks["soma"] > 0], masks["soma"][masks["soma"] > 0]):
        raise RuntimeError("Refinement changed a protected original label")
    return refined, report


def classify_cells(frame):
    result = frame.copy()
    result["morphotype"] = "unclassified_model_unsupported"
    result["mitral_probability"] = np.nan
    values = np.column_stack(
        [
            np.log(result.soma_area_um2),
            np.log(pd.to_numeric(result.primary_nucleus_area_um2, errors="coerce")),
            np.log1p(result.cytoplasmic_contrast),
            np.log(result.primary_DAPI_median.clip(lower=1e-9)),
        ]
    )
    eligible = (
        result.analysis_ROI_fraction.to_numpy() >= 0.1
    ) & ~result.touches_frame.to_numpy()
    complete = (
        np.isfinite(values).all(1) & eligible & result.measurement_valid.to_numpy()
    )
    model, info = fit_model(values[complete], FEATURES)
    if model is None or not info["supported"]:
        return result, info
    result.loc[eligible, "morphotype"] = "unclassified_missing_features"
    probabilities, positive = predict(values[complete], model)
    result.loc[complete, "mitral_probability"] = probabilities
    result.loc[complete, "morphotype"] = np.where(
        positive,
        "mitral_like",
        np.where(probabilities <= 0.2, "other_morphotype", "uncertain"),
    )
    marginal = (
        eligible
        & ~complete
        & np.isfinite(values[:, [0, 2]]).all(1)
        & result.measurement_valid.to_numpy()
    )
    for index in np.flatnonzero(marginal):
        row = result.iloc[index]
        p = marginal_membership(row.soma_area_um2, row.cytoplasmic_contrast, model)
        result.loc[result.index[index], "mitral_probability"] = p
        result.loc[result.index[index], "morphotype"] = (
            "mitral_like"
            if p >= 0.8 and row.soma_area_um2 >= model[-1]
            else "other_morphotype" if p <= 0.2 else "uncertain"
        )
    return result, info
