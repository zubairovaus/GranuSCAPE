# Release checklist (maintainer action required)

This archive is a cleaned **release candidate**, not yet a formally licensed or externally validated release. The Python algorithms have not been modified during the packaging cleanup.

- [ ] Select and add a public `LICENSE` file after confirming rights and agreement of code contributors.
- [ ] Confirm author attribution, contacts, affiliation and citation format.
- [ ] Confirm redistribution rights for the four bundled `.ilp` trained ilastik projects (including any underlying training images/labels embedded in them).
- [ ] Run tests in the intended Windows environment: `python -m pip install -e ".[test]"` and `python -m pytest -q`.
- [ ] Complete one end-to-end sample analysis using supported CZI, Cellpose and ilastik environments; check output QC.
- [ ] Confirm correct dependencies, pinned environment versions and model identifiers for manuscript reproducibility.
- [ ] Confirm Monte Carlo settings and results against the archived manuscript analysis.
- [ ] Review included `examples/selected_cells_136.csv` for permission to share study metadata.
- [ ] Check file contents and repository history for confidential material before making repository public.
- [ ] Set repository visibility public, tag an immutable release, archive through Zenodo and cite resulting DOI.

No license is included pending the maintainer's decision. Until a license is granted, the public repository will generally be viewable but does not automatically confer permission to reuse or redistribute its contents.
