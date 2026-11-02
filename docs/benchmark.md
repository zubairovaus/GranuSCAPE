# Runtime benchmark protocol

The benchmark performs fresh, complete CZI processing. Existing segmentations and previous HDF5 results are not reused. Each run writes a new output folder. The benchmark is separate from read-only QC regeneration.

## Recommended sequence

1. First verify the benchmark command on the existing 3–5-image test subset. For the manuscript, freeze approximately 10 representative images spanning sparse/dense fields and different soma/punctum counts. Select them before examining runtime; retain failures and do not select only the fastest images.
2. Use the same workstation, GPU/CPU mode and final analysis configuration throughout. Connect a laptop to AC power, prevent sleep and close unrelated compute jobs. Record the storage medium, power mode and relevant conditions manually. The tool records hardware where available but does not control power or thermal throttling.
3. Preinstall dependencies and make the requested Cellpose model available before the benchmark. Downloads and software installation are not the intended workload. Review the first run if an unexpected model download occurs.
4. In test.yaml, point input to the fixed test-image folder and use the intended production settings. For the current default output this means A–D plus cell_ids.png. Set keep_work: false if that is the intended production policy; changing cleanup or PNG settings changes the workload and must be reported.
5. Run BENCHMARK.bat, or use:

```powershell
python -m granuscape benchmark --config test.yaml --repeats 3 --warmup 1
```

The output is a new timestamped benchmark folder next to the configured output directory. `--output` may specify another new/empty folder. For N input images, the default performs 3N + 1 complete image processings. The additional first-image run is retained separately as warm-up. Original CZI files and previous results are preserved. Output HDF5 and QC files require additional disk space.

## Timing scope

Monotonic elapsed time is measured with perf_counter. The main processing stages include native CZI import, input hashing, inference, downstream analysis, HDF5 export, QC PNG export, CSV checkpoints and configured work cleanup. Each image is processed sequentially.

Cellpose and ilastik start fresh subprocesses for every image. Their measured wall times include process startup, model/project loading, execution and completion; they are not pure GPU kernel times. Parent-package import, initial hardware inventory and initial manifest hashing occur before measured run_batch calls. Batch wall time includes input validation and final report refreshes. The per-image elapsed timer ends before the final self-report CSV refresh; that refresh remains included in batch wall time. No physical disk-cache flushing or fsync-to-media benchmark is performed.

The warm-up run conditions filesystem/runtime caches but does not leave the model resident for subsequent subprocesses. Its timing is retained, not silently removed. The three measured repetitions use fixed-seed shuffled image orders. Operating-system caches remain active; describe these as repeated processing measurements rather than forced cold-disk timings.

Detailed segmentation and analysis substeps are nested inside inference_seconds and analysis_seconds. Do not add nested substeps to their parent values. Top-level costs plus other_seconds form a per-image decomposition; batch totals also include batch-level work.

## Outputs

- benchmark_protocol.json: fixed configuration, source-code/model hashes and timing scope.
- benchmark_inputs.csv: image names, original paths, file sizes and hashes.
- environment.json: Python/numerical libraries, OS, logical CPUs, available hardware inventory, requested GPU mode, thread environment and disk space. RAM values are installed capacity, not peak use.
- backend_versions.json: actual backend metadata from a completed HDF5 when available, including Cellpose/Torch/device information and ilastik version where exposed.
- benchmark_images.csv: every image timing, including warm-up, failures, dimensions and object counts.
- benchmark_batches.csv: complete run wall times and successful batch throughput.
- benchmark_per_image.csv: within-image medians over measured repeats, restricted to images that completed every planned repeat.
- benchmark_summary.csv: median, quartiles and range of those image-level medians, with top-level/nested scope labels.
- benchmark_status.json: completeness and failure accounting. primary_summary_eligible is false if any planned image lacks all successful measured repetitions.
- BENCHMARK_REPORT.zip: compact report files, without raw CZI, HDF5, probability arrays or PNGs.

If warm-up fails, the benchmark stops so an invalid setup is not repeatedly timed. In measured runs, failures remain in the report. Do not publish an eligible-only summary as a complete-panel benchmark when primary_summary_eligible is false. Technical completion does not establish correct segmentation; expert quality review remains separate.

## Presentation in the manuscript

Report the number of distinct images and repeats separately, image dimensions/calibration and workload counts, hardware/software, GPU mode, QC outputs and the timing boundaries. Summarize processing time per image by first taking its median across repeats, then reporting the median and IQR across images. Batch throughput is based on measured complete-batch wall time, not 3600 divided by a median image time. A useful figure shows the top-level stage decomposition and individual-image runtime versus soma/punctum counts. Warm-up/first-run timings and failures belong in the accompanying report or supplement.

These measurements characterize the selected workload on the tested workstation. They do not by themselves establish segmentation accuracy, peak RAM/VRAM requirements or runtime for a much larger batch. Consolidated table export can grow with batch size.

The native CZI/Cellpose/ilastik benchmark must be run on the target workstation. Synthetic tests validate timing/accounting and numerical parity; they are not a measured full-pipeline performance result for the manuscript.
