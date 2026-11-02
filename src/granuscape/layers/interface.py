import numpy as np
from scipy.ndimage import label, binary_fill_holes, map_coordinates
from scipy.interpolate import BSpline
from skimage.filters import threshold_otsu
from skimage.measure import find_contours
from skimage.morphology import isotropic_closing


def largest_component(mask):
    labels, count = label(mask)
    if not count:
        raise ValueError("No component in density mask")
    sizes = np.bincount(labels.ravel())
    sizes[0] = 0
    return labels == sizes.argmax()


def choose_sparse_interface(density, px, downsample=4, gap_closing_um=0):
    if not np.isfinite(density).all() or density.max() <= density.min():
        raise ValueError("No usable nuclear-density contrast")
    level = float(threshold_otsu(density))
    small = density[::downsample, ::downsample]
    sparse = largest_component(small < level)
    dense_mask = small >= level
    if gap_closing_um:
        radius = gap_closing_um / (px * downsample)
        pad = int(np.ceil(radius)) + 2
        padded = np.pad(dense_mask, pad, mode="edge")
        dense_mask = isotropic_closing(padded, radius)[pad:-pad, pad:-pad]
    components, component_count = label(dense_mask)
    sizes = np.bincount(components.ravel())
    sizes[0] = 0
    ordered = np.argsort(sizes[1:])[::-1] + 1
    dominant_id = int(sizes.argmax())

    def collect(component_id, component_rank):
        region = binary_fill_holes(components == component_id)
        found = []
        for index, yx in enumerate(find_contours(region.astype(np.float32), 0.5)):
            if len(yx) < 30 or np.linalg.norm(yx[0] - yx[-1]) < 2:
                continue
            length = float(np.linalg.norm(np.diff(yx, axis=0), axis=1).sum())
            span = float(np.linalg.norm(yx[0] - yx[-1]))
            if span < 0.25 * min(small.shape):
                continue
            tangent = np.gradient(yx, axis=0)
            normal = np.c_[-tangent[:, 1], tangent[:, 0]]
            normal /= np.maximum(np.linalg.norm(normal, axis=1)[:, None], 1e-08)
            offset = max(2.0, 12.0 / (px * downsample))
            samples = [
                map_coordinates(
                    sparse.astype(float),
                    (yx + sign * offset * normal).T,
                    order=0,
                    mode="constant",
                    cval=0,
                )
                for sign in (-1, 1)
            ]
            adjacency = float(np.maximum(*samples).mean())
            if adjacency < 0.5:
                continue
            score = adjacency * span
            found.append(
                (
                    score,
                    index,
                    yx,
                    adjacency,
                    length,
                    span,
                    component_id,
                    component_rank,
                )
            )
        return found

    candidates = collect(dominant_id, 0)
    used_component_fallback = not bool(candidates)
    if used_component_fallback:
        for rank, component_id in enumerate(ordered):
            if int(component_id) != dominant_id:
                candidates.extend(collect(int(component_id), int(rank)))
    if not candidates:
        raise ValueError(
            "No extended open sparse/dense interface; anatomical review needed"
        )
    candidates.sort(key=lambda row: row[0], reverse=True)
    (
        score,
        index,
        yx,
        adjacency,
        length,
        span,
        selected_component_id,
        selected_component_rank,
    ) = candidates[0]
    dense = binary_fill_holes(components == selected_component_id)
    raw = yx[:, ::-1] * downsample
    h, w = density.shape
    margin = np.minimum.reduce(
        [raw[:, 0], w - 1 - raw[:, 0], raw[:, 1], h - 1 - raw[:, 1]]
    )
    inside = np.flatnonzero(margin > 3 / px)
    if len(inside) < 20:
        raise ValueError("Interface too short inside frame")
    raw = raw[inside[0] : inside[-1] + 1]
    info = dict(
        density_level=level,
        contour_index=index,
        gap_closing_um=gap_closing_um,
        sparse_area_fraction=float(sparse.mean()),
        dense_area_fraction=float(dense.mean()),
        sparse_adjacency_fraction=adjacency,
        candidate_count=len(candidates),
        selection="largest sparse region / dominant dense region open interface",
        geometry_model="quadratic graph in density-interface PCA coordinates",
    )
    if used_component_fallback:
        info.update(
            component_fallback=True,
            selected_dense_component_rank=selected_component_rank,
            selected_dense_component_label=selected_component_id,
            fallback_reason="Dominant dense component has no accepted open interface adjacent to largest sparse region",
            selection="extended open interface of another dense component adjacent to largest sparse region",
            review_required=True,
        )
    return (raw, info)


def fit_open_arc(raw_um):
    raw = np.asarray(raw_um, dtype=float)
    origin = raw.mean(axis=0)
    _, _, vh = np.linalg.svd(raw - origin, full_matrices=False)
    along, across = (vh[0], vh[1])
    s = (raw - origin) @ along
    if s[-1] < s[0]:
        along = -along
        s = -s
    low, high = np.quantile(s, [0.01, 0.99])
    if high - low <= 0:
        raise ValueError("Degenerate interface")
    u = (s - low) / (high - low)
    transverse = (raw - origin) @ across
    edges = np.linspace(0, 1, 81)
    points = []
    for a, b in zip(edges[:-1], edges[1:]):
        keep = (u >= a) & (u < b)
        if keep.any():
            points.append(
                (float(np.median(u[keep])), float(np.median(transverse[keep])))
            )
    if len(points) < 10:
        raise ValueError("Insufficient interface extent")
    uv = np.asarray(points)
    a0, a1, a2 = np.polynomial.polynomial.polyfit(uv[:, 0], uv[:, 1], 2)
    c0 = origin + low * along + a0 * across
    c1 = (high - low) * along + a1 * across
    c2 = a2 * across
    controls = np.array([c0, c0 + c1 / 3, c0 + 2 * c1 / 3 + c2 / 3, c0 + c1 + c2])
    model = BSpline([0, 0, 0, 0, 1, 1, 1, 1], controls, 3)
    fitted = model(np.clip(u, 0, 1))
    residual = np.linalg.norm(fitted - raw, axis=1)
    if np.median(residual) > 10 or np.quantile(residual, 0.95) > 20:
        raise ValueError(
            "Smooth arc deviates too far from density interface; review required"
        )
    return (
        model,
        dict(
            interface_fit_median_distance_um=float(np.median(residual)),
            interface_fit_p95_distance_um=float(np.quantile(residual, 0.95)),
        ),
    )


def validate_open_curve(curve):
    p = np.asarray(curve, float)
    if len(p) < 3 or not np.isfinite(p).all():
        raise ValueError("Invalid MCL geometry")
    if np.linalg.norm(p[0] - p[-1]) < 2:
        raise ValueError("MCL must be an open curve")
    segment = np.diff(p, axis=0)
    if np.any(np.linalg.norm(segment, axis=1) < 1e-09):
        raise ValueError("Degenerate MCL segment")

    def cross(a, b):
        return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]

    for i in range(len(segment) - 2):
        a, b = (p[i], p[i + 1])
        c, d = (p[i + 2 : -1], p[i + 3 :])
        o1, o2 = (cross(b - a, c - a), cross(b - a, d - a))
        o3, o4 = (cross(d - c, a - c), cross(d - c, b - c))
        box = (
            np.maximum(np.minimum(a, b), np.minimum(c, d))
            <= np.minimum(np.maximum(a, b), np.maximum(c, d)) + 1e-08
        ).all(axis=1)
        if np.any(box & (o1 * o2 <= 0) & (o3 * o4 <= 0)):
            raise ValueError("MCL self-intersection detected")
    return dict(open_curve=True, self_intersections=0)
