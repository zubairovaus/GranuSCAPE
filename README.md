
<p align="center">
  <img src="branding/granuscape-logo.png"
       alt="GranuSCAPE — Spatial Cellular Analysis and Puncta Exploration"
       width="700">
</p>

# GranuScape

**Layer-aware cellular morphometry and spatial profiling of fluorescent granules.**

GranuScape processes calibrated three-channel olfactory bulb CZI images using Cellpose, ilastik pixel classification, anatomical layer detection, soma refinement and spatial measurements. The package distinguishes marker-positive puncta from a biological claim that every punctum is a stress granule.

## First Windows test

See `docs/windows_test.md`. `SETUP_TEST.ps1` creates an isolated environment and `test.yaml`; `CHECK_TEST.bat` validates settings and `RUN_TEST.bat` runs only the CZI files placed in `data/test`.

## Installation

Use Python 3.11 or 3.12 in a virtual environment:

```powershell
python -m pip install -e .
```

Cellpose and ilastik run in separate existing environments. Their executable paths are supplied in the configuration. The Cellpose model identifier defaults to `cpsam_v2`, matching the source workflow; its availability must be verified in the configured Cellpose environment. Cellpose weights are not bundled. Four existing ilastik projects are included unchanged: two pixel classifiers and two threshold configurations. The threshold projects do **not** perform object classification.

Copy `examples/experiment.yaml`, set the input/output folders and executable paths, and optionally provide an image metadata CSV. Relative paths are resolved against the configuration file, not the working directory. Input folders are searched recursively. Duplicate image stems are rejected.

```powershell
granuscape check --config experiment.yaml
granuscape run --config experiment.yaml
```

The check command validates configuration structure and executable paths; it does not certify backend compatibility. Use a new empty output folder for each run. Source images are read only. Failures are recorded per image and produce a nonzero final exit code; other images continue.

## Modules

| Module | Responsibility |
|---|---|
| `io` | Configuration, native CZI channels, acquisition XML and HDF5 tables |
| `segmentation` | Isolated Cellpose inference, ilastik probabilities and native threshold operator |
| `layers` | Nuclear-density interface, constrained curve and adaptive MCL region |
| `cells` | Expert-authorized merges, signal-based contour expansion, nucleus association and morphotype candidates |
| `granules` | Separate marker objects, soma ownership and inter-marker overlap |
| `measurements` | Cell coverage, perinuclear/distal coverage and spatial descriptors |
| `qc` | Four configurable QC views, numbered overlays and read-only HDF5 re-export |
| `pipeline` | Batch orchestration, failure isolation, provenance and result export |

## Inputs

Only single-scene, single-Z, single-T, native uint16 CZI images with isotropic XY pixels are currently supported. No implicit Z projection or channel identity guessing is performed. The configured CZI channel names must occur exactly once.

The current biological mapping is AF488 → eIF3η, Rhoda → G3BP, DAPI → nuclei. These identities must match the actual staining and acquisition. Changing a name in YAML does not adapt the trained pixel classifiers to a different assay.

Optional metadata columns are `image_id,group,animal,side,slice,grade`. `image_id` equals the CZI filename without its extension. `animal_uid` combines group and animal; animal numbers are not assumed to be globally unique. Grade is descriptive metadata, not a hidden selection filter. Every CZI in the selected folder is processed.

## Outputs

```text
results/
  hdf5/<image_id>.h5
  qc/<image_id>/A_overlay.png
  qc/<image_id>/B_masks.png
  qc/<image_id>/C_cells.png
  qc/<image_id>/D_puncta.png
  qc/<image_id>/cell_ids.png
  cells.csv
  puncta_g3bp.csv
  puncta_eif3.csv
  images.csv
  run_summary.csv
  configuration.json
```

The four default views are A (all contours over locally enhanced RGB), B (translucent masks on black), C (soma/nucleus contours over RGB) and D (puncta and pixel overlap on black). A separate `cell_ids.png` retains the original soma IDs in the MCL region. All somata have one inspection color; color does not imply morphotype or acceptance. Both MCL curves use the same dusty rose color, with a solid interface and dashed ROI borders at 3 native pixels.

Choose views using the optional `qc` section in `examples/experiment.yaml`. `qc.include_cell_ids` defaults to true, independently of the selected views. Existing configurations without a `qc` section receive these defaults. Each QC folder includes `display_settings.json`. See `docs/qc.md` for display settings and re-export from completed HDF5, without CZI import or inference.

Intermediate backend arrays are retained inside HDF5. Temporary work is removed after successful publication unless `keep_work: true`; failed work and logs are retained. Temporary backend storage is redirected to the output folder to reduce use of the Windows system drive.

Each HDF5 stores raw planes, unmodified acquisition XML, original and refined segmentation, MCL geometry, probability arrays, threshold masks, tables and provenance. It does not duplicate the original CZI bytes or trained ilastik project bytes. File/model hashes and exact numerical dependency versions are retained. Acquisition exposure and time fields remain available in the original XML; the core does not apply cross-image exposure normalization.

Every extant soma appears in `cells.csv`, including outside-region and flagged profiles. The quantitative MCL subset requires `analysis_ROI_fraction >= 0.1`, `touches_frame == False`, `measurement_valid == True`, plus visual acceptance of the image. Subset `morphotype == mitral_like` only when the local classification is supported. A complete file is **not** automatically an expert-accepted result.

Object tables retain unassigned and outside-MCL objects for audit. Regional object summaries count objects with any positive MCL intersection; regional areas are clipped to the MCL mask. Cell object ownership requires a unique majority of at least half the object pixels. Coverage uses all positive pixels clipped to the soma, independently of ownership. All pixel area fractions use the whole projected soma as denominator. Nuclear-projected signal is retained and measured separately.

## Review and scope

Results have status `complete_requires_review`. Review the MCL and numbered soma overlays before quantitative use. Nuclear-dominated masks are flagged conservatively when nuclear area occupies at least 80% of the soma. This new diagnostic is a review trigger, **not a validated quality classifier**. Undefined spatial measures are missing, never filled with zero; zero-puncta cells remain in the tables.

Every existing soma seed is evaluated for refinement, regardless of size, morphotype, preliminary contrast or MCL membership. Frame-truncated cells remain excluded. There is no 45 µm² area cutoff or 20-pixel cytoplasm gate in contour refinement. A nonempty extranuclear seed sample and an adequate local background are still necessary; the signal sample size is exported. The area criterion for anatomical soma support is a separate MCL-fitting rule, not a contour-refinement gate.

Signal-based expansion cannot recover a boundary without sufficient extranuclear signal. A skipped expansion is reported explicitly in the `expansion_status` and `expansion_error` columns. `soma_seed_area_um2`, `expansion_applied`, `expansion_signal_pixels` and `expansion_review_flags` support direct auditing of small cells. Unstable proposals are rejected, retaining the original seed. Generic automatic merging and recovery of Cellpose-missed cells are not implemented in this candidate. Expert merges and explicit MCL boundaries can be supplied using the image-hash-bound format in `docs/expert_edits.md`.

Local morphotype classification requires sufficient complete profiles and a supported two-component mixture. Unsupported images receive `unclassified_model_unsupported` rather than an apparently valid zero-mitral result. There is no silently substituted reference model. Anatomical mitral identity is not established by morphology alone.

Group statistics, PCA/UMAP, manual CellCounter validation and manuscript-specific layouts belong in separate analysis workflows. Monte Carlo spatial null models are implemented in the optional `granuscape spatial` stage and are not executed by the core batch command.

## Runtime measurement

`BENCHMARK.bat` repeats fresh complete CZI processing from test.yaml three times and records one separate first-image warm-up run. It writes a new benchmark folder and a compact BENCHMARK_REPORT.zip. No previous HDF5 is reused or overwritten. See docs/benchmark.md for image selection, hardware conditions, timing boundaries and reporting of repeats/failures.

```powershell
python -m granuscape benchmark --config test.yaml --repeats 3 --warmup 1
```

Regular run_summary.csv now includes monotonic stage timings and image workload sizes. Segmentation/analysis component timers are nested; do not add them to their parent totals. Read-only QC-only times remain a separate measurement.

## Validation

```powershell
python -m pip install -e ".[test]"
python -m pytest -q
```

See `docs/validation.md` for performed checks and remaining release requirements. This is an integration candidate, not a scientifically validated final release. Author attribution, redistribution permissions for trained projects and the public code license must be settled before publication.


## Projection-aware spatial analysis

New strict/extended nuclear-distance measurements and optional geometric Monte Carlo models are described in [SPATIAL_ANALYSIS_METHODS_UPDATE.md](SPATIAL_ANALYSIS_METHODS_UPDATE.md). Use RUN_SPATIAL.bat for read-only reanalysis of final HDF5 files. Old spatial fields are preserved and have different primary-only exclusion semantics.

## Public release status

This repository snapshot is version `0.1.0rc5`. For licensing, trained-model redistribution and release validation requirements, see [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md). The bundled example configuration contains workstation-specific example backend paths and must be adapted before running.
