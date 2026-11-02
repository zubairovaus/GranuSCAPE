# Quality-control images

By default each image generates four views plus a separate numbered overlay:

| File | Purpose |
|---|---|
| `A_overlay.png` | All segmentation contours over dimmed, locally enhanced RGB |
| `B_masks.png` | Translucent somata/nuclei and marker masks on black |
| `C_cells.png` | Soma and nucleus boundaries over RGB |
| `D_puncta.png` | Separate marker masks and exact positive-pixel intersection on black |
| `cell_ids.png` | A-style overlay with original soma IDs on dark label backgrounds |

The ID overlay includes every existing soma with analysis_ROI_fraction >= 0.1, irrespective of morphotype or measurement validity. It does not renumber cells. Labels are for inspection; a labeled cell is not automatically accepted for quantitative analysis. If IDs are not needed, set include_cell_ids to false. This switch is independent of views; views: [] with include_cell_ids: true generates only the numbered PNG.

## Configuration

```yaml
qc:
  views: [A, B, C, D]
  include_cell_ids: true
  normalization: clahe
  background_brightness: 0.38
  clahe_kernel_px: [128, 128]
  clahe_clip_limit: 0.005
  clahe_gamma: 1.3
  mcl_color: '#C86F80'
  mcl_line_width_px: 3
```

To select fewer views, use e.g. views: [A, B]. Set normalization to percentile to omit CLAHE and gamma. There is one shared MCL color and width setting: the central interface is solid; ROI borders are dashed. Soma contours are gold, nuclei light blue, G3BP coral, eIF3η green and exact punctum overlap white.

Display normalization is fitted on complete channels before composition, not separately on crops. Percentile clipping and CLAHE can suppress low display values or enhance noise. They do not change raw arrays, object masks, classification or any quantitative metric. The black-background views do not require fluorescence normalization. One-pixel emphasis on soma/punctum edges is display only. Scale bars use the original pixel calibration and are placed at bottom right. Legends sit above the native image field. Fonts use DejaVu Sans or Arial where available; an installed fallback may have different glyph metrics.

## Re-export without inference

```powershell
python -m granuscape qc --h5 "path\to\final.h5" --output "path\to\new_qc_folder"
python -m granuscape qc --h5 "path\to\final_hdf5_folder" --output "path\to\another_new_qc_folder" --config "examples\qc.yaml"
```

The qc command opens final GranuScape HDF5 files read-only. It does not require a CZI folder, Cellpose, ilastik or backend executable paths. The optional config may be a full experiment file or a small YAML containing only qc options. Output folders must be new/empty to avoid mixing old and new views. On Windows, RENDER_QC.bat prompts for input and output folders.

The root qc_report.csv records per-file completion/failure. Each image folder includes display_settings.json and, for re-export, source.json with the source HDF5 hash. Existing HDF5 and CSV results are not modified. If a source is incomplete/incompatible, its failure is recorded rather than silently omitted.
