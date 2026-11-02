"""Command line entry point."""

import argparse
import json
from granuscape import __version__
from granuscape.io.config import load_config


def main():
    parser = argparse.ArgumentParser(
        prog="granuscape",
        description="Layer-aware cellular morphometry and granule profiling",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    for command in ("check", "run"):
        subparser = subparsers.add_parser(command)
        subparser.add_argument("--config", required=True)
    qc = subparsers.add_parser(
        "qc", help="Regenerate QC from existing HDF5, without inference"
    )
    qc.add_argument(
        "--h5",
        required=True,
        help="One final HDF5 file or a folder of final HDF5 files",
    )
    qc.add_argument("--output", required=True, help="New, empty QC output folder")
    qc.add_argument(
        "--config",
        help="Optional YAML containing qc settings; backend paths are not required",
    )
    spatial = subparsers.add_parser(
        "spatial",
        help="Projection-aware measurements from final HDF5 without inference",
    )
    spatial.add_argument("--h5", required=True)
    spatial.add_argument("--output", required=True)
    spatial.add_argument(
        "--config", help="YAML with spatial options only, no backend paths"
    )
    spatial.add_argument("--selection", help="CSV with image_id or Snap and soma_id")
    bench = subparsers.add_parser(
        "benchmark", help="Repeat complete CZI processing with stage timings"
    )
    bench.add_argument("--config", required=True)
    bench.add_argument(
        "--output", help="New empty benchmark folder; timestamped by default"
    )
    bench.add_argument("--repeats", type=int, default=3)
    bench.add_argument("--warmup", type=int, default=1)
    bench.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    try:
        if args.command == "spatial":
            import yaml
            from pathlib import Path
            from granuscape.spatial.export import reanalyze

            cfg = (
                yaml.safe_load(Path(args.config).read_text(encoding="utf-8"))
                if args.config
                else {}
            )
            return reanalyze(
                args.h5, args.output, (cfg or {}).get("spatial", cfg), args.selection
            )
        if args.command == "qc":
            from granuscape.qc.export import export_qc, load_display_config

            return export_qc(args.h5, args.output, load_display_config(args.config))
        config = load_config(args.config)
        if args.command == "benchmark":
            from granuscape.benchmark.run import benchmark

            return benchmark(config, args.output, args.repeats, args.warmup, args.seed)
        if args.command == "check":
            print(json.dumps(config, indent=2))
            print(
                "Configuration and executable paths are valid. Inference compatibility still requires a test image."
            )
            return 0
        from granuscape.pipeline.run import run_batch

        return run_batch(config)
    except Exception as error:
        parser.exit(1, f"GranuScape: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
