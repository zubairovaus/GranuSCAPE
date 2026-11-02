import numpy as np
from scipy import ndimage as ndi


def measure_profile(s, n, af, dapi, px, sid):
    mask = s == sid
    ys, xs = np.where(mask)
    pad = round(12 / px)
    box = np.s_[
        max(0, ys.min() - pad) : min(s.shape[0], ys.max() + pad + 1),
        max(0, xs.min() - pad) : min(s.shape[1], xs.max() + pad + 1),
    ]
    local = mask[box]
    dist = ndi.distance_transform_edt(~local) * px
    bg = (dist >= 3) & (dist <= 12) & (s[box] == 0) & (n[box] == 0)
    cyto = local & (n[box] == 0)
    z = np.nan
    if bg.sum() >= 30 and cyto.sum() >= 20:
        v = af[box]
        b = np.median(v[bg])
        mad = np.median(abs(v[bg].astype(float) - b))
        z = max(0, (np.median(v[cyto]) - b) / max(1.4826 * mad, 0.01 * abs(b), 1e-09))
    depth = ndi.distance_transform_edt(np.pad(local, 1))[1:-1, 1:-1] * px
    cy, cx = ndi.center_of_mass(mask)
    radius = np.sqrt(mask.sum() * px * px / np.pi)
    ranks = []
    for nid in np.unique(n[mask]):
        if not nid:
            continue
        nm = n == nid
        inside = nm & mask
        ny, nx = ndi.center_of_mass(nm)
        coverage = inside.sum() / nm.sum()
        d = float(np.median(depth[(n[box] == nid) & local]))
        central = np.exp(-((np.hypot(nx - cx, ny - cy) * px / radius) ** 2))
        size = min(nm.sum() / (0.25 * mask.sum()), 1)
        score = (
            0.4 * coverage
            + 0.3 * d / max(depth.max(), 1e-09)
            + 0.2 * central
            + 0.1 * size
        )
        ranks.append(
            dict(
                nucleus_id=int(nid),
                coverage=float(coverage),
                area_um2=float(nm.sum() * px**2),
                score=float(score),
                median_depth_um=d,
            )
        )
    ranks.sort(key=lambda x: x["score"], reverse=True)
    valid = [r for r in ranks if r["coverage"] >= 0.5 and r["area_um2"] >= 10]
    top = valid[0] if valid else None
    margin = (
        valid[0]["score"] - valid[1]["score"] if len(valid) > 1 else 1 if valid else 0
    )
    accepted = bool(top and top["score"] >= 0.55 and (margin >= 0.15))
    from skimage.measure import regionprops

    props = regionprops(mask.astype("uint8"))[0]
    features = dict(
        soma_solidity=float(props.solidity),
        soma_eccentricity=float(props.eccentricity),
        soma_equiv_diameter_px=float(props.equivalent_diameter_area),
        cytoplasm_area_um2=float(np.sum(mask & (n == 0)) * px**2),
        soma_area_um2=float(mask.sum() * px**2),
        soma_area_px=int(mask.sum()),
        soma_centroid_x_px=float(cx),
        soma_centroid_y_px=float(cy),
        cytoplasmic_contrast=float(z),
        primary_nucleus_id=int(top["nucleus_id"]) if accepted else "",
        proposed_nucleus_id=int(top["nucleus_id"]) if top else "",
        primary_nucleus_area_um2=top["area_um2"] if accepted else "",
        primary_nucleus_accepted=str(accepted),
        primary_nucleus_status=(
            "recomputed_geometry_primary"
            if accepted
            else "recomputed_geometry_ambiguous"
        ),
        primary_DAPI_median=(
            float(np.median(dapi[n == top["nucleus_id"]])) if accepted else np.nan
        ),
        primary_nucleus_score_margin=float(margin),
        primary_nucleus_score=top["score"] if top else "",
    )
    return (
        features,
        dict(
            soma_id=sid,
            ranked_nuclei=ranks,
            accepted_primary=accepted,
            score_margin=float(margin),
            score_weights=dict(
                coverage=0.4, relative_depth=0.3, centrality=0.2, size=0.1
            ),
            score_threshold=0.55,
            margin_threshold=0.15,
            AF488_or_DAPI_intensity_used_for_ranking=False,
        ),
    )
