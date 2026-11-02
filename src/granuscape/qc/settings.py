"""Validated display options, independent of segmentation and measurement settings."""

from copy import deepcopy
import math
import re

DEFAULTS = {
    "views": ["A", "B", "C", "D"],
    "include_cell_ids": True,
    "normalization": "clahe",
    "background_brightness": 0.38,
    "clahe_kernel_px": [128, 128],
    "clahe_clip_limit": 0.005,
    "clahe_gamma": 1.3,
    "mcl_color": "#C86F80",
    "mcl_line_width_px": 3,
}


def resolve_options(options=None):
    if options is None:
        options = {}
    if not isinstance(options, dict):
        raise ValueError("qc must be a mapping")
    unknown = set(options) - set(DEFAULTS)
    if unknown:
        raise ValueError(f"Unknown qc options: {sorted(unknown)}")
    result = deepcopy(DEFAULTS)
    result.update(deepcopy(options))
    views = result["views"]
    if not isinstance(views, list) or any(
        not isinstance(v, str) or v not in "ABCD" or len(v) != 1 for v in views
    ):
        raise ValueError("qc.views must be a list containing A, B, C or D")
    if len(set(views)) != len(views):
        raise ValueError("qc.views contains duplicates")
    if not isinstance(result["include_cell_ids"], bool):
        raise ValueError("qc.include_cell_ids must be true or false")
    if result["normalization"] not in ("clahe", "percentile"):
        raise ValueError("qc.normalization must be clahe or percentile")
    for key, lower, upper in (
        ("background_brightness", 0, 1),
        ("clahe_clip_limit", 0, 1),
        ("clahe_gamma", 0, 10),
    ):
        value = result[key]
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not lower < value <= upper
        ):
            raise ValueError(
                f"qc.{key} must be finite, greater than {lower} and no greater than {upper}"
            )
    kernel = result["clahe_kernel_px"]
    if (
        not isinstance(kernel, (list, tuple))
        or len(kernel) != 2
        or any(type(v) is not int or v < 2 for v in kernel)
    ):
        raise ValueError("qc.clahe_kernel_px requires two integers >= 2")
    result["clahe_kernel_px"] = list(kernel)
    if not isinstance(result["mcl_color"], str) or not re.fullmatch(
        r"#[0-9a-fA-F]{6}", result["mcl_color"]
    ):
        raise ValueError("qc.mcl_color must be a six-digit hexadecimal color")
    width = result["mcl_line_width_px"]
    if type(width) is not int or width < 1 or width > 15 or width % 2 != 1:
        raise ValueError("qc.mcl_line_width_px must be an odd integer from 1 to 15")
    return result
