# HDF5 structure

| Path | Content |
|---|---|
| `/raw/dapi`, `/raw/af488`, `/raw/rhoda` | Native fluorescence planes |
| `/acquisition/metadata_xml` | Original CZI acquisition XML |
| `/acquisition` attribute `pixel_size_um` | Validated XY calibration |
| `/segmentation/nuclei`, `/segmentation/soma` | Unmodified Cellpose label arrays |
| `/segmentation/g3bp`, `/segmentation/eif3` | Marker object labels |
| `/segmentation/*_probabilities` | Full pixel probability arrays |
| `/segmentation/*_high`, `/segmentation/*_low` | Native threshold operator outputs |
| `/cells/labels_merged` | Labels after explicit expert merges |
| `/cells/labels_refined` | Final protected/refined soma labels |
| `/mcl/boundary_xy_px` | Layer interface in native x,y pixels |
| `/mcl/region`, `/mcl/density_core` | Analysis region and density-only band |
| `/mcl/band_parameters` | Arc coordinate, center x/y, normal x/y, core sparse/dense widths and ROI sparse/dense widths |
| `/tables/cells` | One column dataset per cell-table field |
| `/tables/puncta_g3bp`, `/tables/puncta_eif3` | Separate object tables |
| `/tables/image` | Regional and image summary |
| `/qc/diagnostics_json` | Layer, expansion, classification and nuclear association diagnostics |
| `/provenance/json` | Source/model hashes, configuration, expert edits and software versions |

All masks share the native y,x image grid. The boundary is in x,y pixels; normal vectors are dimensionless and widths are in micrometres. Boolean regions are distinct from integer instance labels. Floating missing values denote undefined quantities.


## Projection-aware products (schema 2)

Additional tables: /tables/spatial_cells, /tables/spatial_object_regions, /tables/spatial_radial_profiles, /tables/spatial_null_summary. Parameters: /spatial/parameters_json. Raw null vectors: /spatial/null/soma_ID/marker/variant/mode/metric. Sidecars have schema granuscape_spatial and source checksum metadata. Legacy metrics retain their definitions. See SPATIAL_ANALYSIS_METHODS_UPDATE.md.
