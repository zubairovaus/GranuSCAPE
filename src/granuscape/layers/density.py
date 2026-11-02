import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from skimage.measure import regionprops_table
from skimage.draw import line


def props_df(labels, prefix):
    d = pd.DataFrame(
        regionprops_table(
            labels,
            properties=[
                "label",
                "area",
                "centroid",
                "equivalent_diameter_area",
                "eccentricity",
                "solidity",
            ],
        )
    )
    return d.rename(
        columns={
            "label": f"{prefix}_id",
            "area": f"{prefix}_area_px",
            "centroid-0": f"{prefix}_centroid_y_px",
            "centroid-1": f"{prefix}_centroid_x_px",
            "equivalent_diameter_area": f"{prefix}_equiv_diameter_px",
            "eccentricity": f"{prefix}_eccentricity",
            "solidity": f"{prefix}_solidity",
        }
    )


def nucleus_density(nuc_labels, sigma_um, px):
    props = pd.DataFrame(
        regionprops_table(nuc_labels, properties=["label", "centroid"])
    )
    point = np.zeros(nuc_labels.shape, dtype=np.float32)
    yy = np.clip(np.rint(props["centroid-0"]).astype(int), 0, nuc_labels.shape[0] - 1)
    xx = np.clip(np.rint(props["centroid-1"]).astype(int), 0, nuc_labels.shape[1] - 1)
    np.add.at(point, (yy, xx), 1.0)
    return gaussian_filter(point, sigma=sigma_um / px)


def rasterize_curve(shape, curve):
    img = np.zeros(shape, np.uint8)
    pts = np.rint(curve).astype(int)
    for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
        rr, cc = line(y0, x0, y1, x1)
        ok = (rr >= 0) & (rr < shape[0]) & (cc >= 0) & (cc < shape[1])
        img[rr[ok], cc[ok]] = 1
    return img
