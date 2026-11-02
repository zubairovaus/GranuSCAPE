"""Projected soma coverage, nuclear-distance profiles and angular polarization."""

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree


def cell_metrics(cells, soma, nuclei, granules, objects, pixel_size):
    result = cells.copy()
    rows = []
    boxes = ndi.find_objects(soma)
    for row in result.itertuples():
        sid, box = row.soma_id, boxes[row.soma_id - 1]
        mask = soma[box] == sid
        area_pixels = int(mask.sum())
        a, b = granules["g3bp"][box] > 0, granules["eif3"][box] > 0
        primary = (
            int(row.primary_nucleus_id)
            if str(row.primary_nucleus_id) not in ("", "nan")
            else 0
        )
        record = {}
        primary_mask = nuclei[box] == primary if primary else np.zeros(mask.shape, bool)
        # Include full primary nuclear geometry for centroid and nearest-boundary distance.
        if primary:
            nuclear_box = ndi.find_objects((nuclei == primary).astype(np.uint8))[0]
            y0, x0 = min(box[0].start, nuclear_box[0].start), min(
                box[1].start, nuclear_box[1].start
            )
            y1, x1 = max(box[0].stop, nuclear_box[0].stop), max(
                box[1].stop, nuclear_box[1].stop
            )
            full = nuclei[y0:y1, x0:x1] == primary
            distance = ndi.distance_transform_edt(~full) * pixel_size
            distance = distance[
                box[0].start - y0 : box[0].stop - y0,
                box[1].start - x0 : box[1].stop - x0,
            ]
            ny, nx = ndi.center_of_mass(full)
            yy, xx = np.indices(mask.shape)
            angles = np.arctan2(
                yy + box[0].start - (ny + y0), xx + box[1].start - (nx + x0)
            )
            perinuclear = mask & ~primary_mask & (distance <= 2)
            distal = mask & ~primary_mask & (distance > 2)
        for name, positive in (
            ("g3bp", a),
            ("eif3", b),
            ("intersection", a & b),
            ("union", a | b),
        ):
            pixels = int(np.count_nonzero(positive & mask))
            record[name + "_area_um2"] = pixels * pixel_size**2
            record[name + "_coverage"] = pixels / area_pixels
            record[name + "_nuclear_projection_area_um2"] = float(
                np.count_nonzero(positive & mask & (nuclei[box] > 0)) * pixel_size**2
            )
            record[name + "_perinuclear_coverage"] = (
                float(positive[perinuclear].mean())
                if primary and perinuclear.any()
                else np.nan
            )
            record[name + "_distal_coverage"] = (
                float(positive[distal].mean()) if primary and distal.any() else np.nan
            )
            angular = np.full(8, np.nan)
            if primary:
                for sector in range(8):
                    domain = (
                        mask
                        & ~primary_mask
                        & (angles >= -np.pi + sector * np.pi / 4)
                        & (angles < -np.pi + (sector + 1) * np.pi / 4)
                    )
                    if domain.any():
                        angular[sector] = positive[domain].mean()
            record[name + "_angular_coverage_cv"] = (
                float(np.std(angular) / np.mean(angular))
                if np.isfinite(angular).all() and angular.mean() > 0
                else np.nan
            )
            if name in objects:
                assigned = objects[name].loc[objects[name].soma_id == sid]
                record[name + "_count"] = len(assigned)
                record[name + "_mean_object_area_um2"] = float(assigned.area_um2.mean())
                coordinates = (
                    assigned[["centroid_x_px", "centroid_y_px"]].to_numpy(float)
                    * pixel_size
                )
                record[name + "_nearest_neighbor_um"] = (
                    float(cKDTree(coordinates).query(coordinates, k=2)[0][:, 1].mean())
                    if len(coordinates) >= 2
                    else np.nan
                )
                if primary and len(coordinates):
                    dx = coordinates[:, 0] - (nx + x0) * pixel_size
                    dy = coordinates[:, 1] - (ny + y0) * pixel_size
                    radius = np.hypot(dx, dy)
                    valid = radius > 0
                    record[name + "_angular_resultant"] = (
                        float(
                            abs(np.mean((dx[valid] + 1j * dy[valid]) / radius[valid]))
                        )
                        if valid.any()
                        else np.nan
                    )
                else:
                    record[name + "_angular_resultant"] = np.nan
        rows.append(record)
    import pandas as pd

    return pd.concat([result.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
