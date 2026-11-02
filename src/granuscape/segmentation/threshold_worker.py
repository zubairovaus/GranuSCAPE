"""Run the installed ilastik threshold operator, without object classification.

Use ilastik's own python.exe. Input probabilities must be float YXC [0,1].
The .ilp is read only; no project loading, training or GUI is involved.
"""

import argparse
import json
import os
import sys
from pathlib import Path

_dll_handles = []
if os.name == "nt" and hasattr(os, "add_dll_directory"):
    for directory in (
        Path(sys.executable).parent,
        Path(sys.executable).parent / "Library" / "bin",
    ):
        if directory.is_dir():
            _dll_handles.append(os.add_dll_directory(str(directory)))
import h5py
import numpy as np
import vigra
import ilastik
from lazyflow.graph import Graph
from ilastik.applets.thresholdTwoLevels.opThresholdTwoLevels import OpThresholdTwoLevels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", required=True)
    parser.add_argument("--probabilities", required=True)
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()
    with h5py.File(args.project, "r") as project:
        group = project["ThresholdTwoLevels"]
        values = {
            name: group[name][()].item()
            for name in (
                "HighThreshold",
                "LowThreshold",
                "MinSize",
                "MaxSize",
                "Channel",
                "CoreChannel",
                "CurOperator",
            )
        }
        values["SmootherSigma"] = {
            axis: float(group["SmootherSigma/" + axis][()]) for axis in "xyz"
        }
    with h5py.File(args.probabilities, "r") as source:
        array = np.asarray(source["exported_data"], dtype=np.float32)
    if array.ndim != 3 or not np.isfinite(array).all():
        raise ValueError(
            "Expected finite YXC probability array, got " + str(array.shape)
        )
    if array.min() < -1e-06 or array.max() > 1 + 1e-06:
        raise ValueError("Probabilities outside [0,1]")
    print("Parameters:", json.dumps(values), flush=True)
    print("Input YXC:", array.shape, flush=True)
    operator = OpThresholdTwoLevels(graph=Graph())
    try:
        for name, value in values.items():
            getattr(operator, name).setValue(value)
        operator.InputImage.setValue(vigra.taggedView(array, "yxc"))
        for kind, slot in [
            ("final", operator.Output),
            ("high", operator.SmallRegions),
            ("low", operator.BigRegions),
        ]:
            labels = np.squeeze(slot[:].wait())
            if labels.shape != array.shape[:2]:
                raise ValueError("Unexpected output shape: " + str(labels.shape))
            if not np.issubdtype(labels.dtype, np.integer):
                raise ValueError("Expected integer labels")
            output = Path(args.output_prefix + "_" + kind + ".h5")
            temporary = output.with_suffix(".h5.tmp")
            with h5py.File(temporary, "w") as target:
                dataset = target.create_dataset(
                    "exported_data", data=labels.astype(np.uint32), compression="gzip"
                )
                dataset.attrs["axistags"] = vigra.defaultAxistags("yx").toJSON()
                target.attrs["parameters_json"] = json.dumps(values)
                target.attrs["backend"] = "installed ilastik OpThresholdTwoLevels"
                target.attrs["ilastik_version"] = str(
                    getattr(ilastik, "__version__", "unavailable")
                )
                target.attrs["source_project"] = str(Path(args.project).resolve())
                target.attrs["source_probabilities"] = str(
                    Path(args.probabilities).resolve()
                )
            os.replace(temporary, output)
            count = len(np.unique(labels)) - int(np.any(labels == 0))
            print(kind + ": " + str(count) + " labels; " + str(output), flush=True)
    finally:
        operator.cleanUp()


if __name__ == "__main__":
    main()
