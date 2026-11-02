"""Batch orchestration with per-image isolation and atomic HDF5 publication."""

from granuscape.benchmark.timing import Timer
from datetime import datetime, timezone
import importlib.metadata
import json
import hashlib
from pathlib import Path
import shutil
import time
import traceback
import warnings
import h5py
import numpy as np
import pandas as pd
from granuscape import __version__
from granuscape.io.czi import read_czi, sha256
from granuscape.io.hdf5 import write_table, write_json
from granuscape.segmentation.backends import infer
from granuscape.pipeline.analysis import analyze_arrays
from granuscape.qc.plots import save_qc


def code_hash():
    root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    for path in sorted(root.rglob("*.py")):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean_json(value.tolist())
    if isinstance(value, (float, np.floating)) and (not np.isfinite(value)):
        return None
    if isinstance(value, np.generic):
        return value.item()
    return value


def run_batch(config, selected_files=None):
    from granuscape.qc.settings import resolve_options

    config = dict(config, qc=resolve_options(config.get("qc")))
    inputs = (
        sorted(Path(config["input"]).rglob("*"))
        if selected_files is None
        else [Path(p).resolve() for p in selected_files]
    )
    inputs = [p for p in inputs if p.is_file() and p.suffix.lower() == ".czi"]
    if not inputs:
        raise ValueError("No CZI files found")
    stems = [p.stem.casefold() for p in inputs]
    if len(stems) != len(set(stems)):
        raise ValueError(
            "Duplicate image stems in input folder; resolve before processing"
        )
    output = Path(config["output"])
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Output folder is not empty; choose a new output folder")
    output.mkdir(parents=True, exist_ok=True)
    metadata = None
    if config.get("metadata_csv"):
        metadata = pd.read_csv(config["metadata_csv"], dtype=str).fillna("")
        required = {"image_id", "group", "animal", "side"}
        if not required <= set(metadata):
            raise ValueError("Metadata CSV requires image_id, group, animal, side")
        if metadata.image_id.duplicated().any():
            raise ValueError("Metadata CSV has duplicate image_id values")
        missing = set((p.stem for p in inputs)) - set(metadata.image_id)
        if missing:
            raise ValueError(f"Missing image metadata: {sorted(missing)}")
        metadata = metadata.set_index("image_id")
    (output / "configuration.json").write_text(
        json.dumps(config, indent=2), encoding="utf-8"
    )
    edits = {}
    if config.get("edits_json"):
        entries = json.loads(Path(config["edits_json"]).read_text(encoding="utf-8"))[
            "images"
        ]
        for entry in entries:
            if entry["image_id"] in edits:
                raise ValueError("Duplicate image_id in expert edits")
            edits[entry["image_id"]] = entry
    images, cells, objects, report = ([], [], {"g3bp": [], "eif3": []}, [])
    spatial_frames = {
        k: [] for k in ("cells", "object_regions", "null_summary", "radial_profiles")
    }
    for index, source in enumerate(inputs, 1):
        clock = time.perf_counter()
        timer = Timer()
        image_id = source.stem
        print(f"{index}/{len(inputs)} {image_id}", flush=True)
        work = output / "work" / image_id
        work.mkdir(parents=True)
        status = {"image_id": image_id, "status": "failed"}
        try:
            attributes = dict(config.get("experiment", {}))
            if metadata is not None:
                attributes.update(metadata.loc[image_id].to_dict())
            attributes["image_id"] = image_id
            if attributes.get("animal") and attributes.get("group"):
                attributes["animal_uid"] = (
                    attributes["group"] + ":" + attributes["animal"]
                )
            with timer.span("input_hash"):
                source_digest = sha256(source)
            edit = edits.get(image_id, {})
            if edit and edit.get("source_sha256") != source_digest:
                raise ValueError("Expert edits are bound to a different CZI hash")
            with timer.span("import"):
                raw, pixel_size, xml = read_czi(source, config["channels"])
            with timer.span("inference"):
                masks, backends = infer(raw, config, work, timer=timer)
            with timer.span("analysis"):
                with warnings.catch_warnings(record=True) as captured:
                    warnings.simplefilter("always")
                    result = analyze_arrays(
                        raw,
                        masks,
                        pixel_size,
                        config["density_sigma_um"],
                        boundary=(
                            np.asarray(edit["mcl_boundary_xy_px"], float)
                            if "mcl_boundary_xy_px" in edit
                            else None
                        ),
                        edits=edit.get("merges", []),
                        timer=timer,
                        spatial=config.get("spatial"),
                        image_id=image_id,
                    )
            result["diagnostics"]["warnings"] = [str(w.message) for w in captured]
            for frame in [
                result["cells"],
                *result["objects"].values(),
                result["spatial"]["cells"],
                result["spatial"]["object_regions"],
                result["spatial"]["null_summary"],
                result["spatial"]["radial_profiles"],
            ]:
                for key, value in attributes.items():
                    frame[key] = value
            image = dict(attributes, **result["image"])
            status.update(
                width_px=raw["dapi"].shape[1],
                height_px=raw["dapi"].shape[0],
                megapixels=raw["dapi"].size / 1000000.0,
                source_bytes=source.stat().st_size,
                soma_count=len(result["cells"]),
                nucleus_count=len(np.unique(masks["nuclei"]))
                - int(np.any(masks["nuclei"] == 0)),
                g3bp_count=len(result["objects"]["g3bp"]),
                eif3_count=len(result["objects"]["eif3"]),
            )
            target = output / "hdf5" / (image_id + ".h5")
            target.parent.mkdir(exist_ok=True)
            temporary = target.with_suffix(".partial.h5")
            provenance = dict(
                package="GranuScape",
                package_version=__version__,
                source_code_sha256=code_hash(),
                created_utc=datetime.now(timezone.utc).isoformat(),
                source_sha256=source_digest,
                expert_edits=edit,
                source_name=source.name,
                configuration=config,
                backends=backends,
                dependencies={
                    name: importlib.metadata.version(name)
                    for name in (
                        "numpy",
                        "scipy",
                        "pandas",
                        "h5py",
                        "scikit-image",
                        "scikit-learn",
                        "czifile",
                    )
                },
            )
            with timer.span("hdf5_export"):
                with h5py.File(temporary, "w") as h5:
                    h5.attrs.update(
                        schema="granuscape",
                        schema_version="2",
                        status="complete_requires_review",
                    )
                    h5.create_dataset(
                        "acquisition/metadata_xml",
                        data=xml,
                        dtype=h5py.string_dtype("utf-8"),
                    )
                    h5["acquisition"].attrs["pixel_size_um"] = pixel_size
                    write_json(h5, "provenance/json", clean_json(provenance))
                    write_json(
                        h5, "qc/diagnostics_json", clean_json(result["diagnostics"])
                    )
                    for name, array in raw.items():
                        h5.create_dataset(
                            "raw/" + name, data=array, compression="gzip", shuffle=True
                        )
                    for name, array in masks.items():
                        h5.create_dataset(
                            "segmentation/" + name,
                            data=array,
                            compression="gzip",
                            shuffle=True,
                        )
                    for name, array in (
                        ("cells/labels_merged", result["labels_merged"]),
                        ("cells/labels_refined", result["labels_refined"]),
                        ("mcl/boundary_xy_px", result["boundary"]),
                        ("mcl/region", result["region"]),
                        ("mcl/density_core", result["core"]),
                        ("mcl/band_parameters", result["band"]),
                    ):
                        h5.create_dataset(
                            name, data=array, compression="gzip", shuffle=True
                        )
                    from granuscape.spatial.export import write_spatial

                    write_spatial(h5, result["spatial"])
                    write_table(h5, "tables/cells", result["cells"])
                    for marker in objects:
                        write_table(
                            h5, "tables/puncta_" + marker, result["objects"][marker]
                        )
                    write_table(h5, "tables/image", pd.DataFrame([image]))
            with timer.span("qc_export"):
                save_qc(
                    raw,
                    result["labels_refined"],
                    masks["nuclei"],
                    masks,
                    result["boundary"],
                    result["region"],
                    result["cells"],
                    pixel_size,
                    output / "qc" / image_id,
                    options=config.get("qc"),
                )
            if result["spatial"]["options"]["diagnostics"]:
                from granuscape.spatial.export import diagnostic

                diagnostic(result["spatial"], output / "qc" / image_id / "spatial")
            for key in ("cells", "object_regions", "null_summary", "radial_profiles"):
                spatial_frames[key].append(result["spatial"][key])
                frame = pd.concat(spatial_frames[key], ignore_index=True)
                if len(frame):
                    frame.to_csv(output / ("spatial_" + key + ".csv"), index=False)
            temporary.replace(target)
            images.append(image)
            cells.append(result["cells"])
            for marker in objects:
                objects[marker].append(result["objects"][marker])
            status.update(status="complete_requires_review", hdf5=str(target))
            if not config.get("keep_work", False):
                with timer.span("work_cleanup"):
                    shutil.rmtree(work)
        except Exception as error:
            status["error"] = str(error)
            (work / "error.log").write_text(traceback.format_exc(), encoding="utf-8")
        report.append(status)
        with timer.span("csv_export"):
            if images:
                pd.DataFrame(images).to_csv(output / "images.csv", index=False)
                pd.concat(cells, ignore_index=True).to_csv(
                    output / "cells.csv", index=False
                )
                for marker, frames in objects.items():
                    pd.concat(frames, ignore_index=True).to_csv(
                        output / ("puncta_" + marker + ".csv"), index=False
                    )
        status.update(
            {name + "_seconds": value for name, value in timer.seconds.items()}
        )
        status["elapsed_seconds"] = time.perf_counter() - clock
        parents = (
            "input_hash",
            "import",
            "inference",
            "analysis",
            "hdf5_export",
            "qc_export",
            "work_cleanup",
            "csv_export",
        )
        status["other_seconds"] = max(
            0.0,
            status["elapsed_seconds"]
            - sum((timer.seconds.get(name, 0.0) for name in parents)),
        )
        pd.DataFrame(report).to_csv(output / "run_summary.csv", index=False)
    failures = sum((r["status"] == "failed" for r in report))
    print(
        f"Finished: {len(report) - failures} completed, {failures} failed. Output: {output}"
    )
    return 1 if failures else 0
