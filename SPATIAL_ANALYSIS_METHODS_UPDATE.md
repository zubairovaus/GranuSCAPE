# Projection-aware nucleus-centered spatial analysis

Implementation: GranuScape 0.1.0rc5; spatial schema 2. All terminology below concerns two-dimensional projections. A nuclear-overlapping punctum is not thereby an intranuclear granule. The supplied study description specifies physical 4 µm widefield tissue sections; that acquisition information is author-provided, not inferred from masks.

## 1. Audited dependency order and compatibility

The uploaded src.zip and the archived full 0.1.0rc4 package contained identical Python source. Existing order: import/inference; explicit merges; MCL interface detection; soma refinement; final nuclear association and morphometry; adaptive ROI; conditional morphotype classification; whole-object ownership; legacy spatial metrics. The update retains these dependencies and adds projection-aware analysis after ownership, using final refined soma labels and accepted primary nucleus IDs. Monte Carlo never influences segmentation, nuclear selection or morphotype.

Legacy tables/cells fields, including *_perinuclear_coverage, *_distal_coverage and raw angular descriptors, retain their original definitions and values. They exclude the primary nucleus only. They are **not** the new all-nuclei-excluded metrics. New definitions reside in tables/spatial_cells and separate CSV files. Original punctum tables gain owned-object spatial annotations without adding or splitting their rows. Main pipeline HDF5 root schema_version becomes 2. Read-only spatial reanalysis creates separate sidecars (schema=granuscape_spatial, schema_version=2), never overwriting the original.

## 2. Exact geometry

Let S be the refined soma, N* its accepted primary nuclear profile, M_N the union of all nuclear labels and M_c the positive mask for a marker. All masks share native y,x coordinates. Primary nuclear pixels take priority if supplied Boolean masks overlap. The integer label representation itself cannot retain separate overlapping nuclear instances or recover multiple physical nuclei merged into one DAPI label; no automatic splitting or ground-truth claim is made. Diagnostics explicitly flag this unresolved limitation and disconnected primary labels.

The distance d(x) is the Euclidean distance between native pixel centers and the nearest primary nuclear pixel center, multiplied by XY calibration. This matches the old EDT convention; it is not an interpolated membrane-boundary distance. Full primary nuclear geometry is included even where N* extends outside S.

Disjoint regions:

- R_nuc = S ∩ N*.
- R_foreign = S ∩ (M_N \ N*).
- R_peri = (S \ M_N) ∩ {0 < d ≤ r}, r=2 µm.
- R_dist = (S \ M_N) ∩ {d > r}.

These four masks are disjoint and their union equals S, verified in code. R_assoc = R_nuc ∪ R_peri is a deliberately overlapping reporting union, not a fifth partition element. A_strict=S\M_N; A_extended=A_strict ∪ R_nuc. Each region's native pixel area and marker-positive area are exported in µm². No region is identified with a true three-dimensional cytoplasmic compartment.

Coverage C_R=|M_c∩R|/|R|. Strict ratio E_peri=C_peri/C_dist and strict difference Δ_peri=C_peri−C_dist. Extended ratio E_assoc=C_assoc/C_dist and difference Δ_assoc=C_assoc−C_dist are named nucleus-associated, not purely cytoplasmic. Differences are exported as fractions (multiply by 100 for percentage points). Empty region: coverage NaN. Zero distal coverage: ratios NaN, with QC, even when numerator is positive; defined differences remain available. No infinity or arbitrary pseudocount is introduced. A zero marker mask gives zero coverage in nonempty regions. Missing primary nucleus gives explicit unavailable records, not zero spatial metrics.

G3BP, eIF3η, positive-pixel intersection and union receive separate coverage measurements. Derived-mask component counts are labeled count_unit=connected_component_of_derived_mask; they are not new biological punctum counts. G3BP/eIF3η counts always retain original instance identity.

## 3. Objects and ownership

The original owner policy is unchanged: unique maximal soma intersection and at least 50% of whole-object area; tied maxima are unassigned. Pixel coverage uses all positive pixels clipped to the reporting region, independently of ownership. tables/spatial_object_regions contains one association per object and intersected soma, so rows from different somata are not independent object counts.

For whole object P, f_nuc=|P∩N*∩S|/|P|, with analogous fractions for foreign, peri and distal regions. Their sum equals the object's soma-intersection fraction; it need not equal one for a boundary-crossing object. Whole-object area, soma-intersection area, overlap areas, owner, ambiguity, class and signed centroid distance are retained. outside_primary_area_um2 is whole-object area minus the primary overlap inside S; it includes pixels outside S. The centroid uses full-object geometry and a bilinear signed EDT sample (inside negative); the binary centroid-in-primary flag uses the nearest pixel. Pixels are not split into separate counted objects.

Exclusive association classes, in priority order:

1. Ambiguous/foreign-nucleus-associated if the unique owner is not this soma (including unassigned) or any foreign-nucleus overlap occurs.
2. Nuclear-overlapping if primary overlap is positive.
3. Perinuclear if outside-nuclear peri pixels exceed distal pixels.
4. Distal if distal pixels exceed peri pixels.
5. Equal peri/distal split: ambiguous_zone_tie.

Thus a nuclear-overlapping object can have both nuclear and peri pixel fractions while retaining one class. Every zone exports unique intersecting IDs and a separate unique-owned count. Associated count is the cardinality of the union, not nuclear count plus peri count. Counts for adjacent zones must not be summed as though disjoint objects.

## 4. Conditional Monte Carlo models

Ordinary runs default to monte_carlo=false. When enabled, default B=1000, seed=17, modes=[full_mask,centroid_only]. Strict and extended are always separate variants with separate seeds and arrays. A SHA-256-derived seed includes global seed, image ID, soma ID, nucleus ID, marker, variant and mode; image processing order does not alter results.

The modeled object cohort is **uniquely owned whole objects**. Its union supplies the matching observed null-statistic value, separately recorded as observed in spatial_null_summary. This is not necessarily identical to all-positive-pixel coverage in spatial_cells, which can include unassigned or neighboring objects. The cohort field prevents confusing these estimands.

### full_mask

Original object coordinates, orientation, area and shape are preserved; there is no clipping, rotation, resizing, fragmentation or replacement by a centroid. All valid integer translations fully contained in the allowed domain are enumerated for each object; one is sampled uniformly and independently in each realization. Overlap between simulated objects is allowed by policy. Object count remains fixed; pixel coverage uses the Boolean union and does not double count overlapping pixels. This is an independent-location model, not a hard-core exclusion model. Config max_candidates=200000 bounds the number of anchor positions per object. Exceeding this resource limit or having no feasible translation yields QC, not approximate biased placement.

**Observed-domain compatibility:** if any object in this fixed cohort crosses the selected allowed domain, the full-mask variant is unavailable with observed_object_crosses_allowed_domain. Whole objects crossing local crop geometry likewise yield QC. They are not clipped and the cohort is not silently reduced to manufacture a null result. Consequently strict nulls are often unavailable for cells containing nuclear-overlapping owned objects; extended nulls permit primary overlap and can remain available. Foreign overlap remains prohibited in both models. This is an explicit conservative feasibility policy, not evidence that the observed objects are incorrect. Observed coverage and object records are retained regardless.

The strict variant evaluates perinuclear coverage, perinuclear enrichment ratio, coverage difference and strict angular coverage CV. Extended evaluates primary nuclear-projection coverage, nucleus-associated coverage, associated enrichment ratio/difference and extended angular coverage CV. Intersection/union nulls operate on connected components of their observed derived masks and do not represent a joint generative interaction model for the two biological markers.

### centroid_only

For angular resultant only, the original number of owned centroids is preserved and integer sites are sampled uniformly with replacement from the allowed domain. Coinciding sites are permitted. Observed centroids must lie in the same domain; otherwise QC reports incompatibility. Shapes and coverage are not inferred from this baseline. At least two non-central valid centroids are required. The angular reference is the full primary nuclear centroid. Angular resultant is |mean(exp(iθ))|. Sector coverage CV uses eight angular sectors; every sector must have available area and the mean coverage must be positive. Full-mask simulation, not centroid sampling, normalizes sector coverage CV.

### Statistics

Per metric: observed, requested count, finite successful count, null mean, sample SD, 2.5/50/97.5 percentiles, midpoint-tie empirical percentile, two-sided equal-tail empirical p, z-score and QC. With B_f finite values,

p = min(1, 2 min[(1+#(T_b≥T_obs))/(B_f+1), (1+#(T_b≤T_obs))/(B_f+1)]).

Undefined ratios are excluded only from that ratio's null summary, with incomplete-defined-replicates QC. All original B positions, including NaNs, remain stored in the HDF5 array. Zero null variance leaves z NaN and is explicitly flagged. No valid simulations or undefined observation leaves p and z NaN. These conditional per-cell p-values do not demonstrate transport mechanisms and require a downstream multiplicity/replication policy before biological hypothesis testing.

## 5. New tables, HDF5 and diagnostics

- tables/spatial_cells → spatial_cells.csv: long cell × marker table; region areas, positive areas, coverages, ratios/differences, unique counts, raw domain-specific angular metrics, primary ID, image ID, geometry/radius/schema and QC. Animal metadata propagate when provided; absent animal identity is not guessed.
- tables/spatial_object_regions → spatial_object_regions.csv: object × intersected soma associations and fractional distributions, original object ID, owner and ambiguity.
- tables/spatial_radial_profiles → spatial_radial_profiles.csv: all-nuclei-excluded coverage in 0.5 µm bins, with explicit zone denominators; empty bin coverage NaN. This field is new and must replace the old primary-only radial export for new figures.
- tables/spatial_null_summary → spatial_null_summary.csv: model-specific long null summaries, observed estimand/cohort, seed, statistics, successes and QC.
- spatial/parameters_json: resolved configuration; spatial attributes identify schema and primary priority.
- spatial/null/soma_ID/marker/strict_or_extended/full_mask_or_centroid_only/metric: raw B-vector, compressed. Arrays are not expanded into CSV columns.
- Separate diagnostic PNGs: disjoint regions, extended associated region, observed puncta, one feasible full-mask realization, and null histograms with observed reference lines. Only a configurable number are rendered per image.

Changes: new spatial/geometry.py, config.py, analysis.py, nulls.py, export.py; integrations in pipeline/analysis.py, pipeline/run.py, io/config.py, cli.py; package version and schema; tests/test_projection_spatial.py. Legacy measurements/profiles.py is unchanged.

## 6. Configuration and use

```yaml
spatial:
  perinuclear_radius_um: 2.0
  monte_carlo: false
  simulations: 1000
  seed: 17
  modes: [full_mask, centroid_only]
  max_candidates: 200000
  diagnostics: true
  diagnostic_cells: 3
```

Main pipeline accepts this spatial block. `examples/spatial.yaml` enables Monte Carlo explicitly for read-only reanalysis. On the workstation use RUN_SPATIAL.bat (existing Python path configurable in the launcher), select the measured_01 HDF5 directory and a new output directory. Supply the frozen selection CSV to restrict to selected somata; do not apply measurement_valid filtering. The wrapper sets PYTHONPATH to this package's src so it does not accidentally import an older installed copy. A permanent install can use `python -m pip install -e .` in this package directory.

```powershell
python -m granuscape spatial --h5 "...\measured_01\hdf5" --output "...\spatial_projection" --config examples\spatial.yaml --selection selected_cells.csv
```

No Cellpose/ilastik inference is performed by this command. Selection requires soma_id plus image_id or Snap for multi-image inputs. Missing selected IDs fail explicitly. Sidecars bind to the source HDF5 SHA-256 and original files remain read-only.

## 7. Validation and current completion boundary

See the bundled validation report and exact test log for final counts. Synthetic tests cover partition areas and scale, primary priority, foreign nuclei, partial primary profiles, shared nuclear/peri objects without count inflation, foreign/primary-spanning objects, soma boundaries, empty/zero regions, missing primary, disconnected/merged-profile limitation, fragmented/asymmetric domains, exact shape placement and failure/resource-limit paths, repeatability, finite-sample p-values, ratios, full and centroid variants, native HDF5 round trip, read-only reanalysis and legacy parity when Monte Carlo is toggled. A substituted native-run integration tests the batch HDF5/CSV/QC path with external backends mocked; it is not fresh CZI inference.

Real integration uses historical Snap-54527 masks, not the current 20-image or 136-cell biological cohort. Twelve large profiles with foreign-nucleus overlap were inspected computationally, producing 24 marker comparisons against the unchanged legacy formulas: 12 perinuclear and 9 distal coverages changed. The source file SHA-256 remained unchanged. A separate native-mask crop (soma 483) exercised B=1000 models; see real_mc_summary.csv for finite versus unavailable statistics. This is a technical check, not biological validation; the historical image overlaps ilastik development material and must not be presented as an independent accuracy evaluation.

The latest 20 final HDF5 files are not available in this environment. Their reanalysis, the full current-cohort old/new comparison and a fresh CZI/Cellpose/ilastik run have **not** been performed here. The new command/launcher is provided for that necessary workstation step. Existing manuscript plots remain results of the earlier geometry and must not be relabeled as the new strict/extended results. Do not claim the entire real-data validation stage is complete before these outputs are checked.

## 8. Methods draft

Following soma refinement, the primary DAPI-positive nuclear profile was assigned using soma-guided geometric association. All segmented nuclear profiles were then used to partition each soma projection into primary nuclear projection, foreign nuclear overlap, available perinuclear (0–2 µm from the primary profile) and distal regions. Pixel coverages were calculated separately for G3BP, eIF3η and their positive-pixel intersection and union. Strict perinuclear enrichment compared available extranuclear perinuclear coverage with distal coverage. An extended nucleus-associated measure additionally included the primary nuclear projection; it was not interpreted as exclusively cytoplasmic enrichment. Nuclear-overlapping puncta retained their original identities and fractional region intersections, without being interpreted as intranuclear granules.

For compatible object cohorts, conditional geometric null distributions were generated using 1000 independent realizations of uniformly sampled valid translations of whole observed object masks within the specified strict or extended domain. Object shapes, areas, orientations and counts were preserved; simulated objects could overlap and coverage used their union. Cohorts incompatible with the domain were reported as unavailable rather than clipped. Angular resultant was evaluated against a separate centroid-location baseline, whereas sector coverage CV used full-mask realizations. Random seeds were deterministically derived for each image, soma, marker and model variant. Null summaries included finite-simulation counts, mean, SD, quantiles, empirical two-sided equal-tail p-values and standardized deviations where defined. These models describe conditional departures from geometrically constrained random placement and do not establish a biological localization mechanism.
