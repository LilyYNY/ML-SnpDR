# ML-SnpDR

ML-SnpDR extends [subnetDR](https://github.com/LilyYNY/subnetDR) with
machine-learning module scoring, module-level survival analysis, drug-response
evidence integration, and an explicit candidate-selection gate.

Only modules listed in `selected_modules.tsv` are allowed to enter the original
subnetDR sequence/SMILES, binding-affinity, and perturbation-response steps.

## Repository status

The repository is being implemented stage by stage.

- [x] Repository metadata, configuration contract, module identifier and dry-run entry point
- [ ] Module manifest adapter for ModuleSelection outputs
- [ ] Core34 feature engineering and nested-OOF ML scoring
- [ ] Survival, drug-response and candidate-selection gate
- [ ] Selected-only sequence/SMILES, binding and perturbation scoring
- [ ] Paper regression tests and end-to-end example

The current version is a bootstrap release and does not yet execute the full
scientific pipeline.

## Planned workflow

```text
DEPs -> Network construction -> Module division -> Size prefilter
     -> Module manifest -> Core34 features -> ML OOF ranking
     -> Survival + drug response -> Candidate selection
     -> selected_modules.tsv
     -> Sequence/SMILES -> Binding affinity -> Perturbation score
```

## Configuration-first dry run

```r
library(MLSnpDR)

run <- run_ml_snpdr_pipeline(
  config = "config/paper_luad.yml",
  dry_run = TRUE
)

run$plan
```

The LUAD paper-reproduction configuration currently encodes:

1. module size greater than 9 for the initial prefilter;
2. top 10 nested-OOF target-subtype probabilities;
3. high-score-worse survival direction and log-rank P below 0.05;
4. module size of at least 30 genes;
5. PRISM drug count and drug-response density for final ranking.

## Documentation

- [Architecture and migration plan](docs/architecture.md)
- [Data policy](data-raw/README.md)
- [Model artifact policy](models/README.md)

## License and attribution

ML-SnpDR and its upstream subnetDR-derived code are distributed under the MIT
License. See [LICENSE.md](LICENSE.md) and [NOTICE](NOTICE).
