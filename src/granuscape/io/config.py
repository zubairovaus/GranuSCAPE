"""Read explicit, portable experiment and backend configuration."""

from pathlib import Path
import yaml

REQUIRED = {"input", "output", "channels", "backends"}


def load_config(path):
    path = Path(path).resolve()
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or not REQUIRED <= config.keys():
        raise ValueError("Configuration requires input, output, channels and backends")
    allowed = REQUIRED | {
        "experiment",
        "metadata_csv",
        "gpu",
        "keep_work",
        "density_sigma_um",
        "edits_json",
        "qc",
        "spatial",
    }
    unknown = set(config) - allowed
    if unknown:
        raise ValueError(f"Unknown configuration fields: {sorted(unknown)}")
    if set(config["channels"]) != {"dapi", "af488", "rhoda"}:
        raise ValueError(
            "channels must map dapi, af488 and rhoda to exact CZI channel names"
        )
    if len(set(config["channels"].values())) != 3:
        raise ValueError(
            "Each biological channel must have a distinct CZI channel name"
        )
    for key in ("input", "output", "metadata_csv", "edits_json"):
        if config.get(key):
            p = Path(config[key]).expanduser()
            config[key] = str((path.parent / p).resolve() if not p.is_absolute() else p)
    if not Path(config["input"]).is_dir():
        raise ValueError("Input folder does not exist")
    if Path(config["input"]) == Path(config["output"]):
        raise ValueError("Input and output folders must differ")
    for key in ("cellpose_python", "ilastik", "ilastik_python"):
        if not Path(config["backends"].get(key, "")).is_file():
            raise ValueError(f"Missing backend executable: {key}")
    sigma = float(config.get("density_sigma_um", 12))
    if sigma <= 0:
        raise ValueError("density_sigma_um must be positive")
    config["density_sigma_um"] = sigma
    from granuscape.qc.settings import resolve_options

    config["qc"] = resolve_options(config.get("qc"))
    from granuscape.spatial.config import options

    config["spatial"] = options(config.get("spatial"))
    return config
