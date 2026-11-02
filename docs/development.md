# Development boundaries

The core package owns per-image inference, anatomical regions, morphology, puncta and quality review. It exports self-contained files and explicit measurements. It does not decide experimental group comparisons, remove animals to improve significance or generate manuscript-specific results.

A separate `granuscape-analysis` project can consume the exported HDF5 and CSV schema for animal-level aggregation, left/right comparisons, resampling, paired marker profiles, spatial null simulations and PCA. A manuscript repository can contain figure scripts, fixed panel selections and final captions. Both should pin the core release and model hashes.

Avoid coupling anatomical decisions to marker-positive object counts. Distinguish technical success, measurable profiles and expert-accepted images. New heuristics require independent validation rather than only implementation tests.
