import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import closing, disk


def expand(af, soma, nuc, sid, px, fraction=0.4):
    old = soma == sid
    if not old.any():
        raise RuntimeError("Missing target mask")
    d = ndi.distance_transform_edt(~old) * px
    bg = (d >= 3) & (d <= 12) & (soma == 0) & (nuc == 0)
    cyto = old & (nuc == 0)
    if bg.sum() < 30 or not cyto.any():
        raise RuntimeError(
            "Insufficient background pixels or no extranuclear seed signal"
        )
    b = float(np.median(af[bg]))
    c = float(np.median(af[cyto]))
    mad = float(np.median(abs(af[bg].astype(float) - b)))
    if c <= b:
        raise RuntimeError("No positive diffuse AF488 contrast")
    level = b + fraction * (c - b)
    smooth = ndi.gaussian_filter(af.astype(float), 0.6 / px)
    domain = (smooth >= level) & (d <= 10)
    domain = closing(domain, disk(max(1, round(0.8 / px)), decomposition="sequence"))
    domain |= old
    domain &= (soma == 0) | old
    lab, _ = ndi.label(domain)
    ids, counts = np.unique(lab[old], return_counts=True)
    target = int(ids[np.argmax(counts)])
    mask = lab == target
    mask = ndi.binary_fill_holes(mask)
    signed = ndi.distance_transform_edt(mask) - ndi.distance_transform_edt(~mask)
    mask = ndi.gaussian_filter(signed, 0.6 / px) > 0
    mask |= old
    mask &= (soma == 0) | old
    if ndi.label(mask)[1] != 1:
        raise RuntimeError("Expansion is disconnected")
    return (
        mask,
        dict(
            AF488_background_median=b,
            AF488_background_MAD=mad,
            AF488_cytoplasm_median=c,
            signal_sample_pixels=int(cyto.sum()),
            background_sample_pixels=int(bg.sum()),
            threshold_fraction=fraction,
            AF488_threshold=level,
            maximum_added_distance_um=(
                float(d[mask & ~old].max()) if (mask & ~old).any() else 0.0
            ),
        ),
    )


def rank_nuclei(mask, nuc, px):
    depth = ndi.distance_transform_edt(np.pad(mask, 1))[1:-1, 1:-1] * px
    cy, cx = ndi.center_of_mass(mask)
    radius = np.sqrt(mask.sum() * px**2 / np.pi)
    ranks = []
    for nid in np.unique(nuc[mask]):
        if not nid:
            continue
        nm = nuc == nid
        inside = nm & mask
        ny, nx = ndi.center_of_mass(nm)
        coverage = inside.sum() / nm.sum()
        area = nm.sum() * px**2
        d = float(np.median(depth[inside]))
        score = (
            0.4 * coverage
            + 0.3 * d / max(depth.max(), 1e-09)
            + 0.2 * np.exp(-((np.hypot(nx - cx, ny - cy) * px / radius) ** 2))
            + 0.1 * min(nm.sum() / (0.25 * mask.sum()), 1)
        )
        clipped = bool(
            np.any(nm[0]) or np.any(nm[-1]) or np.any(nm[:, 0]) or np.any(nm[:, -1])
        )
        ranks.append(
            dict(
                nucleus_id=int(nid),
                area_um2=float(area),
                overlap_um2=float(inside.sum() * px**2),
                coverage=float(coverage),
                score=float(score),
                touches_tile_edge=clipped,
            )
        )
    ranks.sort(key=lambda x: x["score"], reverse=True)
    valid = [
        r
        for r in ranks
        if not r["touches_tile_edge"] and r["coverage"] >= 0.5 and (r["area_um2"] >= 10)
    ]
    margin = (
        valid[0]["score"] - valid[1]["score"] if len(valid) > 1 else 1 if valid else 0
    )
    accepted = bool(valid and valid[0]["score"] >= 0.55 and (margin >= 0.15))
    return dict(
        projected_nuclei=ranks,
        proposed_primary_nucleus_id=valid[0]["nucleus_id"] if valid else None,
        accepted_primary=accepted,
        score_margin=float(margin),
        intensities_used_for_ranking=False,
    )
