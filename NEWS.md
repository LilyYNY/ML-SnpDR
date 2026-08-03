# MLSnpDR 0.0.1

- Established the local repository and R/Python package skeleton.
- Added versioned YAML configuration files.
- Added canonical module identifier helpers.
- Added a dry-run stage registry and configuration validation.
- Added a validated ModuleSelection adapter that creates per-module node and
  edge files, SHA256 checksums, `module_manifest.tsv` and group-level QC.
- Added a subnetDR-compatible `module_selection()` that reads ModuleDivision
  directly while preserving the original aggregate outputs.
- Replaced the provisional 15-stage plan with steps 1-6, inserted steps
  6A-6C, and selected-only subnetDR steps 7-9.
- Added explicit per-stage input/output contracts and a configuration-driven
  runner that can execute implemented stages independently.
- Implemented all-module functional annotation and four-panel drug-response
  standardization with coverage audits and canonical DRN handoff files.
- Implemented Core34 feature adaptation, repeated nested-OOF Gradient Boosting
  scoring, per-subtype Top-K selection, and ML gate auditing.
- Implemented survival/drug-evidence triage with a complete stepwise filtering
  trace and one self-contained selected module per subtype.
- Reworked subnetDR steps 7-9 to read only `selected_modules.tsv`, with explicit
  sequence/SMILES, binding-score, sensitivity, and final-candidate contracts.
- Added `run_ML_SnpDR()` as the continuous step-4-through-step-9 entry point.
