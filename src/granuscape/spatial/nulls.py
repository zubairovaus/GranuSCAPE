"""Independent uniform translations of whole masks; overlaps are allowed.

No rotations, clipping, resizing or placement by centroid is used in full_mask.
Every valid integer translation is equiprobable for each object. Inter-object
mask overlap is allowed; coverage uses their union, counts retain object IDs.
"""

import hashlib
import numpy as np
from scipy import ndimage as ndi
from .geometry import field_metrics


def stable_seed(seed, *keys):
    return int.from_bytes(
        hashlib.sha256(("|".join(map(str, (seed,) + keys))).encode()).digest()[:8],
        "little",
    )


def resultant(coords, center):
    if len(coords) < 2:
        return np.nan
    dy, dx = (coords - np.asarray(center)).T
    radius = np.hypot(dx, dy)
    valid = radius > 0
    return (
        float(abs(np.mean((dx[valid] + 1j * dy[valid]) / radius[valid])))
        if valid.sum() >= 2
        else np.nan
    )


def summarize(observed, values, requested):
    values = np.asarray(values, float)
    v = values[np.isfinite(values)]
    out = dict(
        observed=observed,
        requested=requested,
        successful=len(v),
        direction="two_sided_equal_tail",
        null_mean=np.nan,
        null_sd=np.nan,
        null_q025=np.nan,
        null_median=np.nan,
        null_q975=np.nan,
        empirical_percentile=np.nan,
        empirical_p=np.nan,
        z_score=np.nan,
        qc="",
    )
    if not np.isfinite(observed):
        out["qc"] = "observed_undefined"
        return out
    if not len(v):
        out["qc"] = "no_valid_null_values"
        return out
    out.update(
        null_mean=float(v.mean()),
        null_sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan,
        null_q025=float(np.quantile(v, 0.025)),
        null_median=float(np.median(v)),
        null_q975=float(np.quantile(v, 0.975)),
        empirical_percentile=float(
            ((v < observed).sum() + 0.5 * (v == observed).sum()) / len(v)
        ),
    )
    out["empirical_p"] = float(
        min(
            1.0,
            2
            * min(
                ((v >= observed).sum() + 1) / (len(v) + 1),
                ((v <= observed).sum() + 1) / (len(v) + 1),
            ),
        )
    )
    if np.isfinite(out["null_sd"]) and out["null_sd"] > 0:
        out["z_score"] = (observed - out["null_mean"]) / out["null_sd"]
    else:
        out["qc"] = "zero_null_variance"
    if len(v) < requested:
        out["qc"] += ";incomplete_defined_replicates"
    return out


def valid_anchors(domain, coords, limit):
    """Exact integer anchors, chunked to avoid a candidate-by-object allocation."""
    shape = coords.max(axis=0) + 1
    if np.any(shape > domain.shape):
        return np.empty((0, 2), int), "object_larger_than_domain"
    # Bounding-box origin need not itself be an occupied mask pixel.
    anchor_shape = np.asarray(domain.shape) - shape + 1
    if int(np.prod(anchor_shape)) > limit:
        return np.empty((0, 2), int), "placement_candidate_limit"
    yy, xx = np.indices(tuple(anchor_shape))
    candidates = np.column_stack([yy.ravel(), xx.ravel()])
    valid = np.ones(len(candidates), bool)
    for y, x in coords:
        valid &= domain[candidates[:, 0] + y, candidates[:, 1] + x]
    return candidates[valid], "" if valid.any() else "no_feasible_translation"


def simulate(
    objects, observed_positive, reg, sectors, center, opts, seed, variant, mode
):
    domain = reg[variant]
    B = opts["simulations"]
    rng = np.random.default_rng(seed)
    flags = []
    observed = field_metrics(observed_positive, reg, sectors)
    coords = np.array([o["coords"].mean(axis=0) for o in objects]).reshape(-1, 2)
    observed["angular_resultant"] = resultant(coords, center)
    names = (
        ["angular_resultant"]
        if mode == "centroid_only"
        else (
            [
                "perinuclear_coverage",
                "perinuclear_enrichment",
                "perinuclear_coverage_difference",
                "strict_angular_coverage_cv",
            ]
            if variant == "strict"
            else [
                "nuclear_coverage",
                "associated_coverage",
                "nucleus_associated_enrichment",
                "nucleus_associated_coverage_difference",
                "extended_angular_coverage_cv",
            ]
        )
    )
    arrays = {n: np.full(B, np.nan) for n in names}
    first = None
    if not objects:
        flags.append("no_objects")
    if not domain.any():
        flags.append("empty_allowed_domain")
    # A conditional fixed-shape model must match its observed cohort exactly.
    if any(not o.get("valid_source", True) for o in objects):
        flags.append("whole_object_leaves_local_geometry")
    if mode == "full_mask" and any(
        not domain[o["coords"][:, 0], o["coords"][:, 1]].all() for o in objects
    ):
        flags.append("observed_object_crosses_allowed_domain")
    if mode == "centroid_only" and len(coords):
        rounded = np.rint(coords).astype(int)
        if not domain[rounded[:, 0], rounded[:, 1]].all():
            flags.append("observed_centroid_outside_allowed_domain")
    if mode == "centroid_only" and len(objects) < 2:
        flags.append("insufficient_objects")
    if not flags:
        placements = []
        if mode == "full_mask":
            for o in objects:
                offsets = o["coords"] - o["coords"].min(axis=0)
                anchors, error = valid_anchors(domain, offsets, opts["max_candidates"])
                if error:
                    flags.append(error)
                    break
                placements.append((offsets, anchors))
        if not flags:
            sites = np.argwhere(domain)
            for b in range(B):
                if mode == "centroid_only":
                    value = {
                        "angular_resultant": resultant(
                            sites[rng.integers(len(sites), size=len(objects))], center
                        )
                    }
                else:
                    positive = np.zeros(domain.shape, bool)
                    for offsets, anchors in placements:
                        xy = offsets + anchors[rng.integers(len(anchors))]
                        positive[xy[:, 0], xy[:, 1]] = True
                    value = field_metrics(positive, reg, sectors)
                    if b == 0:
                        first = positive.copy()
                for n in names:
                    arrays[n][b] = value[n]
    stats = []
    for name, arr in arrays.items():
        stats.append(
            dict(
                metric=name,
                **summarize(observed[name], arr, B),
                model_qc=";".join(flags),
                mode=mode,
                variant=variant,
                seed=str(seed),
                objects=len(objects),
                successful_configurations=B if not flags else 0,
            )
        )
    return stats, arrays, first
