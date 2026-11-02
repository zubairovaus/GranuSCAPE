#!/usr/bin/env python3
"""Run inside a separate Cellpose environment; isolated Cellpose inference."""

import hashlib
import argparse, json, importlib.metadata
from pathlib import Path
import h5py, numpy as np


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-h5", required=True)
    p.add_argument("--output-h5", required=True)
    p.add_argument("--model", default="cpsam_v2")
    p.add_argument("--gpu", action="store_true")
    args = p.parse_args()
    from cellpose import models
    import torch

    if args.gpu and not torch.cuda.is_available():
        raise RuntimeError(
            "--gpu requested but CUDA is unavailable; check the Cellpose environment or use CPU"
        )
    with h5py.File(args.input_h5, "r") as h:
        dapi = h["raw/dapi"][:]
        af = h["raw/af488"][:]
    zero = np.zeros_like(dapi)
    inputs = {
        "nuclei": np.stack([dapi, zero, zero], axis=-1),
        "soma": np.stack([af, dapi, zero], axis=-1),
    }
    model = models.CellposeModel(gpu=args.gpu, pretrained_model=args.model)
    weights = getattr(model, "pretrained_model", None)
    weights_hash = None
    if isinstance(weights, (str, Path)) and Path(weights).is_file():
        digest = hashlib.sha256()
        with Path(weights).open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        weights_hash = digest.hexdigest()
    results = {}
    for name, a in inputs.items():
        print("Cellpose inference:", name, flush=True)
        masks, _, _ = model.eval(
            a,
            channel_axis=2,
            normalize={"percentile": [1, 99]},
            cellprob_threshold=0.0,
            flow_threshold=0.4,
        )
        results[name] = np.asarray(masks, dtype=np.uint32)
    with h5py.File(args.output_h5, "w") as h:
        h.attrs["weights_sha256"] = weights_hash or "unresolved_requires_verification"
        h.attrs["cellpose_version"] = importlib.metadata.version("cellpose")
        h.attrs["torch_version"] = str(torch.__version__)
        h.attrs["torch_cuda_version"] = str(torch.version.cuda)
        h.attrs["device"] = torch.cuda.get_device_name(0) if args.gpu else "CPU"
        h.attrs["settings_json"] = json.dumps(
            dict(
                model=args.model,
                gpu=args.gpu,
                normalize_percentile=[1, 99],
                flow_threshold=0.4,
                cellprob_threshold=0,
                nuclei_channels=["DAPI", "zero", "zero"],
                soma_channels=["AF488", "DAPI", "zero"],
            )
        )
        for name, a in results.items():
            h.create_dataset(name + "/labels", data=a, compression="gzip")
    print("Saved Cellpose masks", args.output_h5, flush=True)


if __name__ == "__main__":
    main()
