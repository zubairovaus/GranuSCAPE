"""Explicit image-bound expert edits; no image-specific rules in the algorithm."""

import numpy as np
from scipy import ndimage as ndi


def apply_merges(labels, edits):
    result = labels.copy()
    for edit in edits:
        source, target = int(edit["source_id"]), int(edit["target_id"])
        if source == target or min(source, target) <= 0:
            raise ValueError("Merge requires two distinct positive labels")
        border = set(
            np.unique(np.r_[result[0], result[-1], result[:, 0], result[:, -1]])
        )
        if source in border or target in border:
            raise ValueError(
                "Frame-truncated cells cannot be merged into the quantitative cohort"
            )
        if not np.any(result == source) or not np.any(result == target):
            raise ValueError("Merge refers to an absent cell label")
        if ndi.label((result == source) | (result == target))[1] != 1:
            raise ValueError("Merge regions must be connected")
        result[result == source] = target
    return result
