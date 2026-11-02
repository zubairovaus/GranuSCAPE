# Validation and release criteria

## Checks performed in this build

All 23 automated tests pass in the recorded build environment. The automated test suite checks module imports; explicit unsupported local classification; object-count versus pixel-area semantics; inclusion of zero-puncta cells; missing spatial values; protected seed labels during refinement; explicit connected expert merges; rejection of frame-truncated merges; self-intersection rejection; channel identity and calibration parsing; HDF5 table round-trip; integrated array analysis with an explicit boundary; and automatic layer geometry on a synthetic nuclear-density interface. Additional checks verify successful expansion of a 9 µm² seed with only eight extranuclear pixels and explicit signal-failure decisions for every seed without a contrast prefilter.

The real-image MCL geometry check is recorded in `real_image_check.json`. The curve was constructed, but its 95th-percentile interface residual exceeded the 20 µm criterion (22.91 µm); this is a flagged review case, not a passing anatomical validation. It uses existing masks and channels. It is not a CZI/backend end-to-end test, a segmentation accuracy benchmark or an estimate of full processing speed. Synthetic data and a single MCL check cannot establish general biological validity.

The batch export test substitutes synthetic import/inference, executes the numerical analysis and PNG/HDF5/CSV exports, and verifies that a failed second image is reported without losing the completed first image. It does not exercise native CZI reading or external inference.

## Before declaring a final release

1. Run fresh native CZI files through the configured Cellpose and ilastik installations on the target workstation. Check channel mapping, calibration, acquisition XML and model hashes.
2. Review difficult low-contrast soma boundaries and failures of local morphotype classification. Quantify contour agreement against independent expert boundaries; do not use successful execution as evidence of segmentation accuracy.
3. Test representative straight, curved, partial and absent MCL fields. Verify that unsuitable frames fail or remain flagged, and that expert boundaries are bound to the correct source hash.
4. Measure wall time separately for import, inference, layer detection, cell refinement, measurements and QC on documented hardware. Record GPU, software versions, image dimensions and object counts. Compare configurations on the same images.
5. Check every exported measurement against small independent mask calculations, and compare CellCounter annotations under a predefined matching and positivity policy.
6. Freeze dependencies and models, select a code license and confirm redistribution rights. Archive source code, configuration, hashes, test report and a small approved example dataset.

The package does not yet provide generic missed-cell recovery, automatic fragment merging or automatic repair of nuclear-dominated soma masks. The conservative 80% nuclear-fraction warning is an unvalidated review heuristic. These limitations are explicit to prevent a technically completed result from being mistaken for validated segmentation.


## QC export validation

Six additional QC tests verify default A–D plus original cell identifiers, selection of views, ID-only output, rejection of invalid display options, display configuration without inference backends, black rendering of constant channels, and read-only HDF5 regeneration with protection of existing output folders. A native-resolution QC regeneration on the supplied Snap-54544 HDF5 produced all five PNGs and passed PNG integrity checks. Its source HDF5 SHA-256 remained unchanged. The measured QC-only time is recorded in qc_real_image_check.json and is not a full pipeline speed benchmark.


## Timing and benchmark validation

Three benchmark tests verify that elapsed times remain recorded for failed stages, repeat aggregation operates within each image and preserves failure accounting, and repetitions use fresh output folders with unchanged input files. The existing batch test verifies top-level and nested elapsed-time fields. A separate numerical parity check compares masks, geometry, cell tables, punctum tables and regional metrics with the previously packaged analysis on the synthetic fixture; results match. See timing_numerical_check.json. Native CZI/backend performance measurements remain to be collected on the target workstation.
