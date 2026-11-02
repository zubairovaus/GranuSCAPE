"""Object ownership and pixel-mask measurements are deliberately independent."""

import numpy as np
import pandas as pd
from skimage.measure import regionprops

OBJECT_COLUMNS = [
    "punctum_id",
    "area_um2",
    "centroid_x_px",
    "centroid_y_px",
    "in_mcl_fraction",
    "soma_id",
    "soma_fraction",
    "touches_frame",
    "raw_mean_intensity",
    "raw_integrated_intensity",
    "overlap_area_um2",
]


def object_table(labels, soma, region, raw, other, pixel_size):
    rows = []
    for prop in regionprops(labels):
        box = prop.slice
        mask = labels[box] == prop.label
        ids, counts = np.unique(soma[box][mask], return_counts=True)
        valid = ids > 0
        owner, fraction = 0, 0.0
        if valid.any():
            positive_ids, positive_counts = ids[valid], counts[valid]
            maximum = positive_counts.max()
            # Equal maxima are explicitly unassigned.
            fraction = float(maximum / prop.area)
            if (positive_counts == maximum).sum() == 1 and fraction >= 0.5:
                owner = int(positive_ids[np.argmax(positive_counts)])
        values = raw[box][mask].astype(float)
        rows.append(
            dict(
                punctum_id=int(prop.label),
                area_um2=float(prop.area * pixel_size**2),
                centroid_x_px=float(prop.centroid[1]),
                centroid_y_px=float(prop.centroid[0]),
                in_mcl_fraction=float(region[box][mask].mean()),
                soma_id=owner,
                soma_fraction=fraction,
                touches_frame=bool(
                    prop.bbox[0] == 0
                    or prop.bbox[1] == 0
                    or prop.bbox[2] == labels.shape[0]
                    or prop.bbox[3] == labels.shape[1]
                ),
                raw_mean_intensity=float(values.mean()),
                raw_integrated_intensity=float(values.sum()),
                overlap_area_um2=float(
                    np.count_nonzero(other[box][mask]) * pixel_size**2
                ),
            )
        )
    return pd.DataFrame(rows, columns=OBJECT_COLUMNS)
