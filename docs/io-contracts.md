# ML-SnpDR Input/Output Contracts

## 1. General rules

- Canonical output tables are UTF-8 TSV files with one header row and case-sensitive column names.
- Inputs support TSV and CSV, plus XLSX where required. Functional gene sets also support GMT.
- Every module-level table must contain a unique `module_uid`: `<network>__<method>__<subtype>__<module>`.
- File paths in manifests use `/`, are relative, and are resolved against the directory containing the manifest.
- Every stage returns an R data frame and stores the primary-output path in an attribute of the returned object.
- Input and output paths are supplied through function arguments or YAML; drive letters are not hard-coded.
- With `strict = TRUE`, identity conflicts, incomplete coverage, non-finite values, and missing handoff files raise errors immediately.
- A new output directory must not already exist, preventing silent overwrites of formal results.

## 2. Interface summary

| Step | Function | Primary input | Primary output |
|---|---|---|---|
| 01 | `run_diff_expr_analysis()` | Expression matrix, subtype table | `differential_expression.tsv` |
| 02 | `run_network_construction()` | Differential expression, PPI index | `network_manifest.tsv` |
| 03 | `module_division()` | Network manifest | `module_division_manifest.tsv` |
| 04 | `module_selection()` | ModuleDivision root | `module_manifest.tsv` |
| 05 | `functional_annotation()` | `module_manifest.tsv`, gene sets | `module_annotation.tsv` |
| 06 | `drug_response_analysis()` | Manifest, drug-sensitivity roots or index | `drug_response_summary.tsv` |
| 06A | `prepare_module_features()` | Manifest, Core34 table | `module_features.tsv` |
| 06B | `run_nested_ml_scoring()` | Features, feature schema | `ml_scores.tsv`, `ml_top10.tsv` |
| 06B | `prepare_ml_scores()` | Features, existing probability table | Same as above |
| 06C | `triage_modules()` | Manifest, ML Top-K, survival, drug summary | `selected_modules.tsv` |
| 07 | `run_SEQCre()` | Selected modules, protein sequences, drug SMILES | `seq_smiles_manifest.tsv` |
| 08 | `predict_BA()` | Selected modules, sequences/SMILES, binding scores or function | `binding_scores.tsv` |
| 09 | `process_prs_dti()` | Selected modules, binding results, sensitivity table or function | `final_candidates.tsv` |

## 3. Step 1: DEPs

### 3.1 Expression input

The first column of the wide table contains unique, non-empty gene/protein identifiers. All remaining columns are numeric samples:

```text
gene  S001  S002  S003
EGFR  10.2  8.7   4.1
TP53  3.2   3.8   5.6
```

When `gene_column` is empty, the first column is used. Supported formats are `.tsv`, `.txt`, `.csv`, and `.xlsx`.

### 3.2 Phenotype input

```text
Sample  Subtype
S001    C1
S002    C1
S003    C2
```

Sample identifiers must be unique and, in strict mode, must all occur in the expression matrix. Every one-versus-rest comparison requires at least two target samples and two remaining samples.

### 3.3 Method

A Wilcoxon rank-sum test and `p.adjust()` are applied separately to each subtype. A gene is tested only when detection rates in both the target and remaining groups meet `detection_threshold`. The `label` rules are:

```text
up:    p_adjust < p_threshold and fold_change > fc_threshold
down:  p_adjust < p_threshold and fold_change < 1 / fc_threshold
other: non_significant
```

### 3.4 `differential_expression.tsv`

```text
subtype, gene, mean_target, mean_other, fold_change, log2_fold_change,
p_value, p_adjust, detection_target, detection_other, label,
n_target, n_other
```

Additional outputs:

- `differential_expression_significant.tsv`
- `differential_expression_summary.tsv`
- `diff_expression_results_all.xlsx` (optional; compatible with subnetDR)

## 4. Step 2: NetworkConstruction

### 4.1 Differential-expression input

The preferred input is the step-1 `differential_expression.tsv`. A subnetDR Excel workbook containing `<subtype>_DiffResults` worksheets is also accepted. By default, only proteins with `label=up` are used; change this behavior with `include_labels`.

### 4.2 PPI index

Minimum format:

```text
network  edge_file
String   PPI/String/string_edges.tsv
```

Complete optional fields:

```text
network, edge_file, node1_column, node2_column,
score_column, score_min, delimiter,
mapping_file, mapping_id_column, mapping_symbol_column, mapping_delimiter
```

- `network` accepts controlled aliases for String, physicalPPIN, and chengF.
- `edge_file` must contain at least two endpoint columns.
- When `score_column` is declared, only rows satisfying `score > score_min` are retained.
- When source endpoints are database identifiers, a mapping file converts them to gene symbols.
- Relative paths are resolved against the directory containing the PPI index.
- `delimiter` and `mapping_delimiter` may be `tab`, `whitespace`, or an explicit delimiter.

Symbol-level edge files can also be supplied directly as a named vector such as `ppi_sources=c(String="...")`.

### 4.3 `network_manifest.tsv`

```text
network_uid, network, subtype, selected_labels, selected_protein_number,
source_edge_number, node_count, edge_count, ppi_file, ppi_sha256,
source_edge_file, analysis_status, status_reason
```

Every `ppi_file` is a canonical `node1,node2` TSV and contains only edges whose endpoints both belong to the selected differential-protein set for that subtype.

## 5. Step 3: ModuleDivision

### 5.1 Input

The preferred input is the step-2 `network_manifest.tsv`, with minimum columns:

```text
network, subtype, ppi_file
```

A legacy `Netconstruct_result/<subtype>/<network>/ppi_<subtype>.txt` root is also supported.

### 5.2 Method

- Louvain: Louvain community detection on an undirected, deduplicated, self-loop-free PPI network.
- WF: modules from edge-betweenness and label propagation are compared with a hypergeometric overlap test, then P values, overlap size, and a fixed seed define non-overlapping consensus modules.

### 5.3 Outputs for each combination

```text
node_Module_<network>_<method>.txt: node,module
edges_<network>_<method>.txt: node1,node2,module
edge_Module_<network>_<method>.txt: node1,node2,module
```

In the complete edge table, `module=0` denotes cross-module edges or edges excluded from WF consensus modules. `edge_Module_*` retains within-module edges only.

### 5.4 `module_division_manifest.tsv`

```text
division_uid, network, method, subtype, module_count, node_count,
edge_count, intramodule_edge_count,
node_module_file, edge_module_file, intramodule_edge_file,
node_sha256, edge_sha256, source_ppi_file,
seed, wf_pvalue_cutoff, analysis_status, status_reason
```

## 6. Step 4: ModuleSelection

### 6.1 Input

`base_input_path` points to the subnetDR step-3 ModuleDivision root. Every subtype/network/method combination requires a node-module table and a module-labeled edge table. Canonical minimum columns are:

```text
node  module
TP53  1
EGFR  2
```

```text
node1  node2  module
TP53   MDM2   1
EGFR   GRB2   2
```

Module 0 denotes cross-module or unassigned edges and does not enter within-module edge files. The filtering rule is `module_size > numberCutoff`.

### 6.2 `module_manifest.tsv`

```text
module_uid, legacy_module_id, network, method, subtype, module,
module_size, edge_count, self_loops_excluded,
node_file, edge_file, prefilter_pass, prefilter_reason,
node_sha256, edge_sha256,
source_module_file, source_node_file, source_edge_file
```

Every `node_file` has `node,module`; every `edge_file` has `node1,node2,module`. Edge files contain within-module edges only and exclude self-loops.

## 7. Step 5: ModuleAnnotation

### 7.1 Custom gene-set input

Minimum long-table columns:

```text
term_id  gene
GO:1     TP53
GO:1     MDM2
```

Optional columns are `database,description`. Each GMT row uses `term_id,description,gene1,gene2,...`. The background-gene file may be a single-column gene table, or a character vector may be passed directly through `background_genes`.

### 7.2 `module_annotation.tsv`

```text
module_uid, network, method, subtype, module,
database, category, term_id, description,
gene_ratio, background_ratio, rich_factor, fold_enrichment,
p_value, p_adjust, q_value, gene_count, gene_ids,
analysis_status, status_reason
```

Refer to the output table for the definitive columns. Modules without significant terms retain a status row so that they are not mistaken for modules that were never analyzed.

Additional outputs:

- `module_annotation_top<N>.tsv`
- `module_annotation_qc.tsv`
- `module_annotation_parameters.tsv`
- Optional per-module annotation files

## 8. Step 6: DrugResponse

### 8.1 Input mode A: panel-named root directories

```r
c(
  PRISM = "path/to/PRISM",
  GDSC1 = "path/to/GDSC1",
  GDSC2 = "path/to/GDSC2",
  CTRP2 = "path/to/CTRP2"
)
```

The directory structure may follow subnetDR's network/method/subtype/module layout. After parsing, it is immediately frozen into an explicit source index.

### 8.2 Input mode B: explicit index

Required columns:

```text
module_uid, drug_panel, drn_file, drn_info_file
```

Optional columns:

```text
drug_level_file, prediction_file
```

Exactly one of the two input modes must be used.

### 8.3 Canonical DRN

The DRN must provide fields that can be normalized to at least:

```text
protein, drug, p_value
```

Optional fields are `p_adjust,effect_size,effect_direction,tested_samples`. DRN-info requires at least:

```text
node, type
```

`type` is either `protein` or `drug`.

### 8.4 `drug_response_summary.tsv`

```text
module_uid, network, method, subtype, module, module_size,
drug_panel, drug_number, tested_drug_number, drug_response_density,
drn_drug_number, drn_edge_number, p_adjust_method, significance_cutoff,
drn_file, drn_info_file, analysis_status, status_reason
```

`drug_response_density = drug_number / module_size`.

Additional outputs:

- `drug_response_hits.tsv`
- `drug_response_source_index.tsv`
- `drug_response_coverage.tsv`
- Standardized per-module DRN/DRN-info files

## 9. Step 6A: ModuleFeatures

### 9.1 Feature input

Each module occupies one row and must include `module_uid`, or module-identity columns that map uniquely to the manifest. Feature columns must contain finite numeric values. Recommended format:

```text
module_uid, module_size, <34 feature columns>
```

If `feature_columns` is not supplied explicitly, the program excludes identity and QC columns, selects the remaining features, and requires their count to equal `expected_feature_count`.

### 9.2 `module_features.tsv`

```text
module_uid, network, method, subtype, module,
<34 features in the order fixed by feature_schema.json>,
feature_missing_count, feature_qc_pass, feature_qc_reason
```

`module_size` can be an explicitly declared Core34 feature; identity columns do not enter the model.

Additional outputs:

- `feature_schema.json`: feature version, count, order, source path, and SHA256.
- `module_features_qc.tsv`: mapping, missingness, and module-size consistency.

## 10. Step 6B: MLScoring

### 10.1 Existing probability-table input

```text
module_uid, prob_C1, prob_C2, prob_C3, prob_C4
```

Identity, predicted-subtype, and model-metadata columns are permitted, but all four probabilities and all manifest modules must be present. Each row's probabilities should sum to 1 within numeric tolerance.

### 10.2 Bundled nested-OOF input

- The step-6A `module_features.tsv`.
- `feature_schema.json` in the same directory.
- Outer/inner CV settings, random seed, Top-K, and `n_jobs`.

`fast` uses the manuscript-selected Gradient Boosting parameters. `paper` runs the complete inner grid within every outer training fold.

### 10.3 `ml_scores.tsv`

```text
module_uid, network, method, subtype, module,
true_subtype, predicted_subtype,
prob_C1, prob_C2, prob_C3, prob_C4,
target_subtype_probability, probability_margin,
rank_in_subtype, score_type, model_version
```

### 10.4 `ml_top10.tsv`

This table contains every column from `ml_scores.tsv` plus:

```text
target_subtype, top_k, ml_gate, ml_gate_reason
```

The actual filename follows `top_k`; for example, `top_k=10` produces `ml_top10.tsv`.

Additional outputs:

- `ml_gate_audit.tsv` or Python fold/metadata files.
- `fitted_model.joblib` when using the bundled Python path.

## 11. Step 6C: ModuleTriage

### 11.1 Survival input

The table contains one row per ML Top-K module. Canonical columns are:

```text
module_uid, n_samples, n_events, matched_gene_fraction,
cox_hr_per_sd, cox_p, cox_fdr, cutpoint,
logrank_p, logrank_fdr, hr_high_vs_low, survival_direction
```

Common aliases are normalized. By default, `survival_direction` must equal `High_score_worse`.

### 11.2 `module_evidence.tsv`

This table joins the ML Top-K, manifest, survival, and primary drug-sensitivity panel one-to-one by `module_uid`, then adds:

```text
prognosis_gate, prognosis_gate_reason,
module_size_gate, module_size_gate_reason,
drug_response_gate, drug_response_gate_reason, all_gates_pass
```

### 11.3 `module_filtering_stepwise.tsv`

```text
module_uid, subtype, gate_order, gate_name,
gate_pass, cumulative_pass, gate_reason
```

### 11.4 `selected_modules.tsv`

```text
module_uid, network, method, subtype, module, module_size,
primary_drug_panel, node_file, edge_file, drn_file, drn_info_file,
target_subtype_probability, probability_margin, ml_rank_in_subtype,
survival_direction, logrank_p, logrank_fdr, hr_high_vs_low,
drug_number, tested_drug_number, drug_response_density, drn_edge_number,
selection_rank, selection_reason
```

By default, the table contains exactly one row per subtype. All four file paths point to self-contained handoff files copied into the step-6C output directory.

## 12. Step 7: SEQCre

### 12.1 Lookup-table input

Protein sequences:

```text
node, sequence
```

Drug structures:

```text
node, smiles
```

Common identifier aliases such as protein/target/gene and drug/compound are recognized.

### 12.2 `seq_smiles_manifest.tsv`

```text
module_uid, subtype, sequence_file, smiles_file,
required_protein_number, matched_protein_number,
required_drug_number, matched_drug_number,
analysis_status, status_reason
```

Each module's `protein_sequences.tsv` contains `node,sequence,sequence_status`; `drug_smiles.tsv` contains `node,SMILES,smiles_status`.

## 13. Step 8: BindingScore

### 13.1 Binding-score input

```text
module_uid, target_name, drug_name, binding_score
```

Every selected-DRN drug-target pair must be unique. Alternatively, pass `predictor(dpi)`, where the DPI contains:

```text
module_uid, target_name, drug_name, target_seq, drug_smiles
```

### 13.2 `binding_scores.tsv`

```text
module_uid, target_name, drug_name, target_seq, drug_smiles,
binding_score, binding_rank, score_direction
```

This stage also writes `binding_score_manifest.tsv` and per-module DPI/binding-score files.

## 14. Step 9: PScore

### 14.1 Target-sensitivity input

```text
module_uid, target_name, sensitivity
```

Alternatively, pass `sensitivity_function(edges, selected_module_row)`. The function returns target and sensitivity columns.

### 14.2 `perturbation_scores.tsv`

This table contains selected-module identity, every binding field from step 8, and:

```text
sensitivity, perturbation_score, perturbation_rank_in_module
```

where:

```text
perturbation_score = binding_score * sensitivity
```

### 14.3 `final_candidates.tsv`

This table contains the top `top_n` rows per subtype after sorting by perturbation score and adds the following column to the perturbation-score table:

```text
final_rank_in_subtype
```

The stage also writes canonical `target_sensitivity.tsv`.

## 15. Top-level runner configuration

`run_ML_SnpDR()` reads external inputs and stage-output directories from `paths` in the configuration. Key external-input fields are:

```text
expression_file
subtype_file
ppi_index_file
differential_expression_dir
network_construction_dir
module_division_dir
annotation_gene_set_file
annotation_background_gene_file
drug_response_input_index_file or drug_response.input_roots
module_feature_input_file
ml_score_input_file (empty runs the bundled nested OOF implementation)
survival_input_file
protein_sequence_file
drug_smiles_file
binding_score_input_file
sensitivity_input_file
```

The defaults are `from="deps"` and `to="perturbation_score"`, so one `run_ML_SnpDR()` call executes steps 1–9 together with the inserted steps 6A–6C. The `differential_expression`, `network_construction`, and `module_division` configuration sections store thresholds and random seeds for steps 1–3. When `from` targets an intermediate stage, the configuration must point to existing upstream primary outputs. `dry_run=TRUE` validates only the configuration and stage interval and does not create scientific results.
