"""Append-only spatial products and read-only reanalysis of published masks."""

from pathlib import Path
import json
import h5py
import numpy as np
import pandas as pd
from granuscape.io.hdf5 import write_table, write_json, read_table
from .analysis import spatial_analysis
from .config import options


def write_spatial(h, result):
    write_json(h, "spatial/parameters_json", result["options"])
    h["spatial"].attrs.update(
        schema_version="2",
        interpretation="two_dimensional_projection",
        primary_priority=True,
    )
    for key in ("cells", "object_regions", "null_summary", "radial_profiles"):
        write_table(h, "tables/spatial_" + key, result[key])
    for path, arrays in result["null_arrays"].items():
        for metric, values in arrays.items():
            h.create_dataset(
                "spatial/null/" + path + "/" + metric,
                data=values,
                compression="gzip",
                shuffle=True,
            )


def diagnostic(result, directory):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    for i, case in enumerate(result["debug"]):
        reg = case["regions"]
        image = np.zeros((*case["soma"].shape, 3))
        palette = {
            "nuclear": (0.45, 0.65, 0.85),
            "foreign": (0.7, 0.45, 0.7),
            "perinuclear": (0.45, 0.75, 0.65),
            "distal": (0.85, 0.8, 0.55),
        }
        for name, color in palette.items():
            image[reg[name]] = color
        fig, ax = plt.subplots(1, 4, figsize=(14, 3), layout="constrained")
        ax[0].imshow(image)
        ax[0].set_title("Primary / foreign / peri / distal")
        ax[1].imshow(reg["associated"], cmap="gray")
        ax[1].set_title("Extended associated zone")
        ax[2].imshow(image * 0.5)
        ax[2].contour(case["observed"], levels=[0.5], colors=["coral"], linewidths=0.7)
        ax[2].set_title("Observed marker projection")
        if case["random"] is not None:
            ax[3].imshow(case["random"], cmap="gray")
            ax[3].set_title("One full-mask null realization")
        else:
            ax[3].imshow(case["primary"], cmap="gray")
            ax[3].set_title("Primary nuclear profile")
        for a in ax:
            a.axis("off")
        fig.suptitle(f'Soma {case["soma_id"]}: {case["marker"]}, {case["variant"]}')
        fig.savefig(directory / f'soma_{case["soma_id"]}_{i}_regions.png', dpi=170)
        plt.close(fig)
        if case["arrays"]:
            fig, axes = plt.subplots(
                1,
                len(case["arrays"]),
                figsize=(4 * len(case["arrays"]), 3),
                squeeze=False,
                layout="constrained",
            )
            for a, (name, values) in zip(axes[0], case["arrays"].items()):
                v = values[np.isfinite(values)]
                a.hist(v, bins=25, color="#378e86")
                a.set_title(name.replace("_", " "), fontsize=8)
                from .geometry import angular_domains, field_metrics

                sectors, _ = angular_domains(case["primary"].shape, case["primary"])
                value = field_metrics(case["observed"], reg, sectors).get(name, np.nan)
                if np.isfinite(value):
                    a.axvline(value, color="coral", label="Observed")
                    a.legend(fontsize=7)
            fig.savefig(directory / f'soma_{case["soma_id"]}_{i}_null.png', dpi=170)
            plt.close(fig)


def reanalyze(source, output, config=None, selection=None):
    from granuscape.io.czi import sha256
    from granuscape import __version__

    source = Path(source)
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Choose a new empty spatial output directory")
    output.mkdir(parents=True, exist_ok=True)
    opts = options(config)
    selected = pd.read_csv(selection) if selection else None
    if selected is not None and not {"soma_id"} <= set(selected):
        raise ValueError(
            "Selection requires soma_id and image_id or Snap for multi-image selection"
        )
    paths = [source] if source.is_file() else sorted(source.rglob("*.h5"))
    if not paths:
        raise ValueError("No HDF5 files")
    if len({p.stem for p in paths}) != len(paths):
        raise ValueError("Duplicate image IDs")
    collected = {
        k: [] for k in ("cells", "object_regions", "null_summary", "radial_profiles")
    }
    audits = []
    matched = set()
    for path in paths:
        image_id = path.stem
        before = sha256(path)
        with h5py.File(path, "r") as h:
            if h.attrs.get("schema") != "granuscape":
                raise ValueError("Requires native GranuScape HDF5 schema")
            cells = read_table(h, "tables/cells")
            if selected is not None:
                if "image_id" in selected:
                    sel = selected[selected.image_id.astype(str) == image_id]
                elif "Snap" in selected:
                    sel = selected[
                        selected.Snap.astype(str) == image_id.replace("Snap-", "")
                    ]
                elif len(paths) == 1:
                    sel = selected
                else:
                    raise ValueError("Multi-image selection requires image_id or Snap")
                if sel.empty:
                    continue
                if not set(sel.soma_id) <= set(cells.soma_id):
                    raise ValueError("Selected soma missing from " + image_id)
                cells = cells[cells.soma_id.isin(sel.soma_id)].copy()
                if sel.soma_id.duplicated().any():
                    raise ValueError("Duplicate selected soma in " + image_id)
                if "primary_nucleus_id" in sel:
                    expected = cells.soma_id.map(
                        sel.set_index("soma_id").primary_nucleus_id
                    )
                    if not np.array_equal(
                        expected.to_numpy(), cells.primary_nucleus_id.to_numpy()
                    ):
                        raise ValueError(
                            "Selected primary nucleus differs from source: " + image_id
                        )
                for column in ("animal", "animal_uid", "group", "side"):
                    if column in sel and column not in cells:
                        cells[column] = cells.soma_id.map(
                            sel.set_index("soma_id")[column]
                        )
                matched.update(sel.index.tolist())
            soma = h["cells/labels_refined"][()]
            nuc = h["segmentation/nuclei"][()]
            granules = {m: h["segmentation/" + m][()] for m in ("g3bp", "eif3")}
            px = float(h["acquisition"].attrs["pixel_size_um"])
            tables = {m: read_table(h, "tables/puncta_" + m) for m in granules}
            result = spatial_analysis(
                cells, soma, nuc, granules, tables, px, opts, image_id
            )
            metadata = {
                k: cells[k].iloc[0]
                for k in ("animal", "animal_uid", "group", "side")
                if k in cells and len(cells)
            }
            for key in collected:
                for k, v in metadata.items():
                    result[key][k] = v
        target = output / (image_id + "_spatial.h5")
        temporary = target.with_suffix(".partial.h5")
        with h5py.File(temporary, "w") as h:
            h.attrs.update(
                schema="granuscape_spatial",
                schema_version="2",
                package_version=__version__,
                source_sha256=before,
                source_name=path.name,
            )
            write_spatial(h, result)
        if opts["diagnostics"]:
            diagnostic(result, output / "qc" / image_id)
        temporary.replace(target)
        if sha256(path) != before:
            raise RuntimeError("Source checksum changed")
        audits.append(
            dict(image_id=image_id, sha256=before, unchanged=True, cells=len(cells))
        )
        for key in collected:
            collected[key].append(result[key])
    if selected is not None and matched != set(selected.index):
        raise ValueError(
            "Some selected rows had no matching source; outputs incomplete"
        )
    for key, frames in collected.items():
        if frames:
            pd.concat(frames, ignore_index=True).to_csv(
                output / ("spatial_" + key + ".csv"), index=False
            )
    pd.DataFrame(audits).to_csv(output / "source_audit.csv", index=False)
    (output / "spatial_configuration.json").write_text(json.dumps(opts, indent=2))
    return 0
