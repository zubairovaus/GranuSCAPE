import numpy as np
from scipy import ndimage as ndi
from .contours import expand


def compute(af, s, n, px, ids, threshold=0.4):
    out = s.astype("uint32", copy=True)
    cost = np.full(s.shape, np.iinfo(np.int32).max, dtype="int32")
    winner = np.zeros(s.shape, "uint32")
    proposal_count = np.zeros(s.shape, "uint16")
    pending = {}
    boxes = ndi.find_objects(s)
    records = []
    pad = int(np.ceil(12 / px))
    shape = s.shape
    for sid in ids:
        if sid <= 0 or sid > len(boxes) or boxes[sid - 1] is None:
            records.append(dict(soma_id=sid, status="absent_label"))
            continue
        b = boxes[sid - 1]
        if (
            b[0].start == 0
            or b[1].start == 0
            or b[0].stop == shape[0]
            or (b[1].stop == shape[1])
        ):
            records.append(dict(soma_id=sid, status="excluded_source_frame_profile"))
            continue
        sl = np.s_[
            max(0, b[0].start - pad) : min(shape[0], b[0].stop + pad),
            max(0, b[1].start - pad) : min(shape[1], b[1].stop + pad),
        ]
        local = s[sl]
        old = local == sid
        if ndi.label(old)[1] != 1:
            records.append(dict(soma_id=sid, status="original_disconnected_review"))
            continue
        try:
            proposal, info = expand(af[sl], local, n[sl], sid, px, threshold)
            strict, _ = expand(af[sl], local, n[sl], sid, px, 0.55)
        except RuntimeError as e:
            records.append(
                dict(soma_id=sid, status="insufficient_signal_review", error=str(e))
            )
            continue
        added = proposal & ~old
        distance2 = np.rint(ndi.distance_transform_edt(~old) ** 2).astype("int32")
        current = cost[sl]
        owner = winner[sl]
        count_local = proposal_count[sl]
        count_local[added] += 1
        ys, xs = np.where(added)
        pending[sid] = (ys + sl[0].start) * shape[1] + xs + sl[1].start
        win = added & (
            (distance2 < current)
            | (distance2 == current) & ((owner == 0) | (sid < owner))
        )
        owner[win] = sid
        current[win] = distance2[win]
        info.update(
            soma_id=sid,
            status="proposed",
            old_area_um2=float(old.sum() * px**2),
            independent_proposed_area_um2=float(proposal.sum() * px**2),
            strict_threshold_IoU=float(
                np.sum(proposal & strict) / np.sum(proposal | strict)
            ),
            overlap_with_other_proposals_pixels=0,
            slice_yxyx=[sl[0].start, sl[0].stop, sl[1].start, sl[1].stop],
        )
        records.append(info)
    out[(s == 0) & (winner > 0)] = winner[(s == 0) & (winner > 0)]
    for r in records:
        if r["status"] != "proposed":
            continue
        sid = r["soma_id"]
        r["overlap_with_other_proposals_pixels"] = int(
            np.sum(proposal_count.ravel()[pending[sid]] > 1)
        )
        y0, y1, x0, x1 = r["slice_yxyx"]
        sl = np.s_[y0:y1, x0:x1]
        local = out[sl]
        original = s[sl] == sid
        lab, _ = ndi.label(local == sid)
        seed = int(np.unique(lab[original])[0])
        islands = (local == sid) & (lab != seed)
        local[islands] = 0
        mask = local == sid
        area = float(mask.sum() * px**2)
        r.update(
            final_area_um2=area,
            area_ratio=area / r["old_area_um2"],
            new_pixels_after_competition=int(np.sum(mask & ~original)),
            removed_disconnected_new_pixels=int(islands.sum()),
            review_flags=[],
        )
        if r["area_ratio"] > 3:
            r["review_flags"].append("large_area_change")
        if r["strict_threshold_IoU"] < 0.75:
            r["review_flags"].append("threshold_sensitive")
        touches = (
            y0 == 0
            and mask[0].any()
            or (x0 == 0 and mask[:, 0].any())
            or (y1 == shape[0] and mask[-1].any())
            or (x1 == shape[1] and mask[:, -1].any())
        )
        if touches:
            r["review_flags"].append(
                "expanded_profile_reaches_frame_exclude_cell_analysis"
            )
        r["status"] = (
            "expanded_requires_review" if mask.sum() > original.sum() else "unchanged"
        )
    if not np.array_equal(out[s > 0], s[s > 0]):
        raise ValueError("Numerical invariant violated")
    return (out, records)
