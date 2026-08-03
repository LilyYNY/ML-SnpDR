# Contributing

## Development rules

1. Keep module-level joins keyed by `module_uid`.
2. Do not add machine-specific absolute paths.
3. Do not mix OOF probabilities with fitted-model probabilities.
4. Preserve the `selected_modules.tsv` boundary for sequence, binding and
   perturbation-response stages.
5. Add or update tests whenever a data contract or selection rule changes.
6. Do not commit patient-level data, API caches or generated result trees.

Run the lightweight checks before proposing changes:

```text
R CMD check --no-manual MLSnpDR_<version>.tar.gz
python -m pytest
```

