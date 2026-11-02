"""Projection-aware spatial measurements, object annotations and null models."""

import numpy as np
import pandas as pd
from scipy import ndimage as ndi
from skimage.measure import regionprops
from .config import options
from .geometry import regions, angular_domains, field_metrics
from .nulls import simulate, stable_seed

MARKERS = ("g3bp", "eif3", "intersection", "union")
REGIONS = ("nuclear", "foreign", "perinuclear", "distal", "associated")


def spatial_analysis(
    cells, soma, nuclei, granules, objects, pixel_size, config=None, image_id=""
):
    opts = options(config)
    radial_rows = []
    rows = []
    associations = []
    null_rows = []
    null_arrays = {}
    debug = []
    positives = {"g3bp": granules["g3bp"] > 0, "eif3": granules["eif3"] > 0}
    positives.update(
        intersection=positives["g3bp"] & positives["eif3"],
        union=positives["g3bp"] | positives["eif3"],
    )
    labels = {
        m: granules[m] if m in ("g3bp", "eif3") else ndi.label(positives[m])[0]
        for m in MARKERS
    }
    # Global objects retain identity and whole-object area, even at soma boundaries.
    props = {m: {p.label: p for p in regionprops(labels[m])} for m in MARKERS}
    owners = {
        m: dict(zip(objects[m].punctum_id, objects[m].soma_id))
        for m in ("g3bp", "eif3")
    }
    from granuscape.granules.objects import object_table

    for m in ("intersection", "union"):
        table = object_table(
            labels[m],
            soma,
            np.ones(soma.shape, bool),
            np.zeros(soma.shape),
            np.zeros(soma.shape, bool),
            pixel_size,
        )
        owners[m] = dict(zip(table.punctum_id, table.soma_id))
    for cell in cells.itertuples():
        sid = int(cell.soma_id)
        nid = int(cell.primary_nucleus_id) if pd.notna(cell.primary_nucleus_id) else 0
        yy, xx = np.nonzero(soma == sid)
        if not len(yy):
            continue
        ny, nx = np.nonzero(nuclei == nid) if nid > 0 else ([], [])
        base = dict(
            image_id=image_id,
            soma_id=sid,
            primary_nucleus_id=nid,
            geometry_model="all_nuclear_profiles_excluded",
            perinuclear_radius_um=opts["perinuclear_radius_um"],
            spatial_schema_version="2",
        )
        if not len(ny):
            for m in MARKERS:
                rows.append(
                    dict(**base, marker=m, spatial_qc="primary_nucleus_missing")
                )
            continue
        y0 = max(0, min(yy.min(), np.min(ny)) - 1)
        y1 = min(soma.shape[0], max(yy.max(), np.max(ny)) + 2)
        x0 = max(0, min(xx.min(), np.min(nx)) - 1)
        x1 = min(soma.shape[1], max(xx.max(), np.max(nx)) + 2)
        # Retain complete intersecting objects, including their external centroids.
        for marker in MARKERS:
            for oid in np.unique(labels[marker][soma == sid]):
                if oid:
                    by0, bx0, by1, bx1 = props[marker][int(oid)].bbox
                    y0, x0 = min(y0, max(0, by0 - 1)), min(x0, max(0, bx0 - 1))
                    y1 = max(y1, min(soma.shape[0], by1 + 1))
                    x1 = max(x1, min(soma.shape[1], bx1 + 1))
        box = np.s_[y0:y1, x0:x1]
        mask = soma[box] == sid
        primary = nuclei[box] == nid
        reg = regions(
            mask, primary, nuclei[box] > 0, pixel_size, opts["perinuclear_radius_um"]
        )
        sectors, center = angular_domains(mask.shape, primary)
        qc = ["two_dimensional_projection", "merged_nuclear_instances_not_resolved"]
        if not np.all(mask[primary]):
            qc.append("primary_nucleus_extends_beyond_soma")
        if reg["foreign"].any():
            qc.append("foreign_nuclear_projection")
        if ndi.label(primary)[1] > 1:
            qc.append("disconnected_primary_nucleus")
        if not reg["perinuclear"].any():
            qc.append("empty_perinuclear_zone")
        if not reg["distal"].any():
            qc.append("empty_distal_zone")
        for marker in MARKERS:
            local = labels[marker][box]
            positive = positives[marker][box]
            rec = dict(
                **base,
                marker=marker,
                **field_metrics(positive, reg, sectors),
                spatial_qc=";".join(qc),
                null_enabled=opts["monte_carlo"],
                count_unit=(
                    "original_marker_object"
                    if marker in ("g3bp", "eif3")
                    else "connected_component_of_derived_mask"
                ),
            )
            domain = reg["strict"]
            radial_index = np.maximum(
                0, np.ceil(reg["distance"] / 0.5 - 1e-10).astype(int) - 1
            )
            for index in range(
                int(radial_index[domain].max()) + 1 if domain.any() else 0
            ):
                zone = domain & (radial_index == index)
                area = int(zone.sum())
                positive_pixels = int((zone & positive).sum())
                radial_rows.append(
                    dict(
                        **base,
                        marker=marker,
                        distance_low_um=index * 0.5,
                        distance_high_um=(index + 1) * 0.5,
                        zone_area_um2=area * pixel_size**2,
                        positive_area_um2=positive_pixels * pixel_size**2,
                        coverage=positive_pixels / area if area else np.nan,
                    )
                )
            for name in REGIONS:
                rec[name + "_area_um2"] = int(reg[name].sum()) * pixel_size**2
                rec[name + "_positive_area_um2"] = (
                    int((reg[name] & positive).sum()) * pixel_size**2
                )
            for name in REGIONS:
                zone_ids = np.unique(local[reg[name]])
                rec[name + "_intersecting_object_count"] = int(
                    np.count_nonzero(zone_ids)
                )
                rec[name + "_owned_object_count"] = sum(
                    owners[marker].get(int(k), 0) == sid for k in zone_ids if k > 0
                )
            if rec["distal_coverage"] == 0:
                rec["spatial_qc"] += ";zero_distal_coverage"
            owned = []
            observed = np.zeros(mask.shape, bool)
            for oid in np.unique(local[mask & positive]):
                prop = props[marker][int(oid)]
                coords = prop.coords - np.array([y0, x0])
                owner = int(owners[marker].get(int(oid), 0))
                inside = mask & (local == oid)
                parts = {
                    n: int(np.count_nonzero(inside & reg[n]))
                    for n in ("nuclear", "foreign", "perinuclear", "distal")
                }
                ambiguous = owner != sid or parts["foreign"] > 0
                if ambiguous:
                    category = "ambiguous_foreign_nucleus_associated"
                elif parts["nuclear"]:
                    category = "nuclear_overlapping"
                elif parts["perinuclear"] > parts["distal"]:
                    category = "perinuclear"
                elif parts["distal"] > parts["perinuclear"]:
                    category = "distal"
                else:
                    category = "ambiguous_zone_tie"
                    ambiguous = True
                cy, cx = np.array(prop.centroid) - np.array([y0, x0])
                iy, ix = np.clip(
                    np.rint([cy, cx]).astype(int), [0, 0], np.array(mask.shape) - 1
                )
                obj = dict(
                    **base,
                    marker=marker,
                    punctum_id=int(oid),
                    owner_soma_id=owner,
                    assignment_ambiguous=ambiguous,
                    spatial_class=category,
                    object_area_um2=prop.area * pixel_size**2,
                    soma_intersection_area_um2=int(inside.sum()) * pixel_size**2,
                    outside_primary_area_um2=(prop.area - parts["nuclear"])
                    * pixel_size**2,
                    centroid_signed_distance_um=float(
                        ndi.map_coordinates(
                            reg["signed_distance"],
                            np.array([[cy], [cx]]),
                            order=1,
                            mode="nearest",
                        )[0]
                    ),
                    centroid_in_primary=bool(primary[iy, ix]),
                    centroid_sampling="bilinear_signed_EDT; nearest-pixel region flag",
                )
                for n, value in parts.items():
                    obj[n + "_overlap_area_um2"] = value * pixel_size**2
                    obj[n + "_fraction_of_whole_object"] = value / prop.area
                obj["associated_fraction_of_whole_object"] = (
                    parts["nuclear"] + parts["perinuclear"]
                ) / prop.area
                associations.append(obj)
                # Include only assigned whole objects; no clipping for shape simulation.
                if owner == sid:
                    fullcoords = prop.coords - np.array([y0, x0])
                    fits = (fullcoords >= 0).all() and (
                        fullcoords < np.array(mask.shape)
                    ).all()
                    owned.append(
                        dict(
                            coords=fullcoords if fits else np.argwhere(inside),
                            valid_source=fits,
                        )
                    )
                    observed[inside] = True
            from .nulls import resultant

            owned_centroids = np.array(
                [o["coords"].mean(axis=0) for o in owned]
            ).reshape(-1, 2)
            for variant in ("strict", "extended"):
                rounded = np.rint(owned_centroids).astype(int)
                valid_centroids = all(o.get("valid_source", True) for o in owned)
                if len(rounded):
                    valid_centroids = valid_centroids and bool(
                        reg[variant][rounded[:, 0], rounded[:, 1]].all()
                    )
                rec[variant + "_owned_angular_resultant"] = (
                    resultant(owned_centroids, center) if valid_centroids else np.nan
                )
                if not valid_centroids:
                    rec["spatial_qc"] += (
                        ";" + variant + "_owned_centroids_outside_domain"
                    )

            rec["owned_object_count"] = len(owned)
            rows.append(rec)
            if opts["monte_carlo"]:
                for variant in ("strict", "extended"):
                    for mode in opts["modes"]:
                        seed = stable_seed(
                            opts["seed"], image_id, sid, nid, marker, variant, mode
                        )
                        stats, arrays, first = simulate(
                            owned,
                            observed,
                            reg,
                            sectors,
                            center,
                            opts,
                            seed,
                            variant,
                            mode,
                        )
                        key = f"soma_{sid}/{marker}/{variant}/{mode}"
                        null_arrays[key] = arrays
                        for stat in stats:
                            null_rows.append(
                                dict(
                                    **base,
                                    marker=marker,
                                    cohort="uniquely_owned_whole_objects",
                                    **stat,
                                )
                            )
                        if first is not None and len(debug) < opts["diagnostic_cells"]:
                            debug.append(
                                dict(
                                    soma_id=sid,
                                    marker=marker,
                                    variant=variant,
                                    mode=mode,
                                    regions=reg,
                                    soma=mask,
                                    primary=primary,
                                    observed=observed,
                                    random=first,
                                    arrays=arrays,
                                )
                            )
        if len(debug) < opts["diagnostic_cells"] and not opts["monte_carlo"]:
            debug.append(
                dict(
                    soma_id=sid,
                    marker="g3bp",
                    variant="observed",
                    mode="none",
                    regions=reg,
                    soma=mask,
                    primary=primary,
                    observed=positives["g3bp"][box],
                    random=None,
                    arrays={},
                )
            )
    # Preserve global object rows; append owned-object spatial fields without duplication.
    assoc = pd.DataFrame(associations)
    annotated = {}
    for m, t in objects.items():
        if len(assoc):
            own = assoc[
                (assoc.marker == m) & (assoc.owner_soma_id == assoc.soma_id)
            ].drop(columns=["image_id", "soma_id", "marker"])
            annotated[m] = t.merge(
                own, on="punctum_id", how="left", validate="one_to_one"
            )
        else:
            annotated[m] = t.copy()
        if "assignment_ambiguous" in annotated[m]:
            annotated[m]["assignment_ambiguous"] = annotated[m][
                "assignment_ambiguous"
            ].map(lambda x: True if pd.isna(x) else bool(x))
        if "spatial_class" in annotated[m]:
            annotated[m]["spatial_class"] = annotated[m]["spatial_class"].fillna(
                "unassigned_no_soma_association"
            )
    return dict(
        cells=pd.DataFrame(rows),
        radial_profiles=pd.DataFrame(radial_rows),
        object_regions=assoc,
        null_summary=pd.DataFrame(null_rows),
        null_arrays=null_arrays,
        objects=annotated,
        debug=debug,
        options=opts,
    )
