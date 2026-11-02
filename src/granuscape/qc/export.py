"""Regenerate QC directly from completed GranuScape HDF5 files, opened read-only."""

import json
from pathlib import Path
import time

import h5py
import pandas as pd
import yaml

from granuscape.io.czi import sha256
from granuscape.io.hdf5 import read_table
from .plots import save_qc
from .settings import resolve_options


def load_display_config(path=None):
    if path is None:
        return resolve_options()
    config = yaml.safe_load(Path(path).read_text(encoding="utf-8-sig"))
    if not isinstance(config, dict):
        raise ValueError("QC configuration must be a mapping")
    return resolve_options(config.get("qc"))


def export_qc(source, output, options=None):
    source, output = Path(source).resolve(), Path(output).resolve()
    options = resolve_options(options)
    files = (
        sorted(
            p for p in source.rglob("*") if p.is_file() and p.suffix.lower() == ".h5"
        )
        if source.is_dir()
        else [source]
    )
    if not files or any(not p.is_file() for p in files):
        raise ValueError("No HDF5 input files found")
    if len({p.stem.casefold() for p in files}) != len(files):
        raise ValueError("Duplicate HDF5 stems; select one final file per image")
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("QC output folder is not empty; choose a new folder")
    output.mkdir(parents=True, exist_ok=True)
    report = []
    for index, path in enumerate(files, 1):
        print(f"QC {index}/{len(files)}: {path.stem}", flush=True)
        clock = time.perf_counter()
        row = {"image_id": path.stem, "status": "failed", "source_h5": str(path)}
        try:
            digest = sha256(path)
            with h5py.File(path, "r") as h5:
                if (
                    h5.attrs.get("schema") != "granuscape"
                    or h5.attrs.get("status") != "complete_requires_review"
                ):
                    raise ValueError("Input is not a completed GranuScape HDF5")
                raw = {key: h5["raw/" + key][:] for key in ("rhoda", "af488", "dapi")}
                soma = h5["cells/labels_refined"][:]
                nuclei = h5["segmentation/nuclei"][:]
                granules = {
                    key: h5["segmentation/" + key][:] for key in ("g3bp", "eif3")
                }
                curve, roi = h5["mcl/boundary_xy_px"][:], h5["mcl/region"][:]
                cells = read_table(h5, "tables/cells")
                pixel_size = float(h5["acquisition"].attrs["pixel_size_um"])
            target = output / path.stem
            outputs = save_qc(
                raw,
                soma,
                nuclei,
                granules,
                curve,
                roi,
                cells,
                pixel_size,
                target,
                options,
            )
            (target / "source.json").write_text(
                json.dumps(
                    {
                        "source_h5": str(path),
                        "source_h5_sha256": digest,
                        "source_modified": False,
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
            row.update(status="complete", png_count=len(outputs))
        except Exception as error:
            row["error"] = str(error)
        row["elapsed_seconds"] = round(time.perf_counter() - clock, 3)
        report.append(row)
        pd.DataFrame(report).to_csv(output / "qc_report.csv", index=False)
    return int(any(row["status"] != "complete" for row in report))
