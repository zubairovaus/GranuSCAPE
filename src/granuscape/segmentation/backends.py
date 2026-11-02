"""External inference in existing, isolated Cellpose and ilastik environments."""

from granuscape.benchmark.timing import Timer
import json
import os
from pathlib import Path
import subprocess
import h5py
import numpy as np
from granuscape.io.czi import sha256

ROOT = Path(__file__).resolve().parent
MODELS = ROOT.parent / "models"


def run_command(command, log, work):
    environment = os.environ.copy()
    temporary = work / "temp"
    temporary.mkdir(exist_ok=True)
    environment.update(TEMP=str(temporary), TMP=str(temporary), TMPDIR=str(temporary))
    with log.open("w", encoding="utf-8") as stream:
        result = subprocess.run(
            list(map(str, command)),
            stdout=stream,
            stderr=subprocess.STDOUT,
            env=environment,
            check=False,
        )
    if result.returncode:
        raise RuntimeError(f"Inference failed; see {log}")


def infer(raw, config, work, timer=None):
    timer = Timer() if timer is None else timer
    work.mkdir(parents=True, exist_ok=True)
    source = work / "raw.h5"
    with h5py.File(source, "w") as h5:
        for name, array in raw.items():
            h5.create_dataset("raw/" + name, data=array)
    backends = config["backends"]
    output = work / "cellpose.h5"
    command = [
        backends["cellpose_python"],
        ROOT / "cellpose_worker.py",
        "--input-h5",
        source,
        "--output-h5",
        output,
        "--model",
        backends.get("cellpose_model", "cpsam_v2"),
    ]
    if config.get("gpu", True):
        command.append("--gpu")
    with timer.span("segmentation_cellpose"):
        run_command(command, work / "cellpose.log", work)
    with h5py.File(output, "r") as h5:
        masks = {key: h5[key + "/labels"][:] for key in ("soma", "nuclei")}
        provenance = {"cellpose": {key: str(value) for key, value in h5.attrs.items()}}
    for marker, channel, classes in (("g3bp", "rhoda", 2), ("eif3", "af488", 3)):
        image = work / (marker + "_raw.h5")
        with h5py.File(image, "w") as h5:
            h5.create_dataset("raw", data=raw[channel])
        probabilities = work / (marker + "_probabilities.h5")
        project = MODELS / (marker + "_pixel.ilp")
        with h5py.File(project, "r") as h5:
            names = [
                v.decode() if isinstance(v, bytes) else str(v)
                for v in h5["PixelClassification/LabelNames"][:]
            ]
        expected = (
            ["punctum", "background"]
            if marker == "g3bp"
            else ["punctum", "diffuse_soma", "background"]
        )
        if names != expected:
            raise ValueError(f"Unexpected model class order: {names}")
        with timer.span("segmentation_" + marker + "_pixel"):
            run_command(
                [
                    backends["ilastik"],
                    "--headless",
                    "--readonly=True",
                    "--project=" + str(project),
                    "--raw_data=" + str(image) + "/raw",
                    "--input_axes=yx",
                    "--export_source=Probabilities",
                    "--output_format=hdf5",
                    "--output_filename_format=" + str(probabilities),
                ],
                work / (marker + "_pixel.log"),
                work,
            )
        with h5py.File(probabilities, "r") as h5:
            array = h5["exported_data"][:]
        if (
            array.shape != raw[channel].shape + (classes,)
            or not np.isfinite(array).all()
        ):
            raise ValueError(
                f"Invalid probability dimensions for {marker}: {array.shape}"
            )
        if (
            array.min() < -1e-06
            or array.max() > 1 + 1e-06
            or (not np.allclose(array.sum(-1), 1, atol=0.02))
        ):
            raise ValueError(f"Invalid probabilities for {marker}")
        threshold = MODELS / (marker + "_threshold.ilp")
        prefix = work / marker
        with timer.span("segmentation_" + marker + "_threshold"):
            run_command(
                [
                    backends["ilastik_python"],
                    ROOT / "threshold_worker.py",
                    "--project",
                    threshold,
                    "--probabilities",
                    probabilities,
                    "--output-prefix",
                    prefix,
                ],
                work / (marker + "_threshold.log"),
                work,
            )
        with h5py.File(str(prefix) + "_final.h5", "r") as h5:
            masks[marker] = h5["exported_data"][:]
            provenance[marker] = {
                "threshold_parameters": json.loads(h5.attrs["parameters_json"]),
                "ilastik_version": str(h5.attrs.get("ilastik_version", "unavailable")),
            }
        masks[marker + "_probabilities"] = array
        for kind in ("high", "low"):
            with h5py.File(str(prefix) + "_" + kind + ".h5", "r") as h5:
                masks[marker + "_" + kind] = h5["exported_data"][:]
    provenance["models_sha256"] = {p.name: sha256(p) for p in MODELS.glob("*.ilp")}
    return (masks, provenance)
