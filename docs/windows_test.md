# First workstation test

1. Extract the archive so that the project root is `D:\sg_local\GranuScape` and contains `pyproject.toml`. Replace the previously extracted source folder with this archive; do not mix files from different builds. Keep the existing analysis and inference environments unchanged.
2. Run in PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File "D:\sg_local\GranuScape\SETUP_TEST.ps1" -Python "D:\sg_local\.venv\Scripts\python.exe"
```

The installer creates a separate `D:\sg_local\granuscape_env`, installs GranuScape and runs its automated checks. It redirects temporary installation files to `D:\sg_local\granuscape_temp` and disables pip's download cache for this installation. Cellpose and ilastik are reused from their existing environments. Internet access is needed to install missing Python dependencies. The PowerShell installer itself has not been executed on Windows in the build environment.

3. Copy 3–5 CZI files into `D:\sg_local\GranuScape\data\test`. Choose one clearly segmented field, one curved MCL field, one low-contrast soma field and, if available, a field with a small soma seed. Do not move the original files or copy the full dataset for the first test.
4. Open `test.yaml`. The automatically written file uses JSON syntax, which is valid YAML. Verify `cellpose_python`, `ilastik` and `ilastik_python`. The default Cellpose path uses your Windows user profile. Check that the `cpsam_v2` model resolves in that environment. Check channel names against the CZI metadata. `gpu` may be changed to false to request CPU inference.
5. Run `CHECK_TEST.bat`; it checks configuration and executable paths. This does not run inference or certify model availability.
6. Run `RUN_TEST.bat`. The first run writes to `D:\sg_local\GranuScape\results\test_01`. For another run, change the output folder in `test.yaml` to a new empty folder, for example `./results/test_02`.
7. Inspect `run_summary.csv` for failures and timing. In `qc\<image_id>`, inspect `A_overlay.png`, `B_masks.png`, `C_cells.png`, `D_puncta.png` and `cell_ids.png`. Check MCL placement, cell borders, protected neighboring cells, puncta and the scale bar. In `cells.csv`, compare `soma_seed_area_um2` with `soma_area_um2`, and inspect `expansion_status`, `expansion_applied`, `expansion_error` and `expansion_review_flags`. A small seed must have an explicit decision, never a size-based skip. Larger growth is not automatically better.
8. Keep `keep_work: true` for this first test so failed inference can be inspected. Share `run_summary.csv`, `images.csv`, `cells.csv`, the used `test.yaml` and the QC PNGs for the test images when reviewing the result. Full HDF5 files are not initially necessary. Backend logs are under `work\<image_id>`.

Every non-frame soma seed is evaluated for refinement. Unchanged masks may be correct. Missing extranuclear signal, insufficient background, disconnected seeds or unstable growth require review. Frame-truncated cells remain excluded from quantitative analysis. Automatic checks do not replace expert inspection or a native CZI/backend validation run.


## Refresh QC from existing results

After replacing the source files, rerun SETUP_TEST.ps1 to install the updated package and its explicit Pillow dependency. Existing test.yaml files remain intact; missing qc settings receive defaults. Do not delete or replace data, results or personal configuration files when extracting the source archive.

Run `RENDER_QC.bat`, provide the existing final HDF5 folder (for example `D:\sg_local\GranuScape\results\test_01\hdf5`) and a new empty output folder (for example `D:\sg_local\GranuScape\results\qc_02`). The launcher uses examples/qc.yaml. Existing HDF5, CSV results, Cellpose and ilastik are not rerun or changed.
