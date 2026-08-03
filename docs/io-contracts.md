# ML-SnpDR 输入输出契约

## 1. 通用规则

- 标准输出表为 UTF-8 TSV，一行表头，列名区分大小写。
- 输入支持 TSV/CSV；需要时支持 XLSX。功能基因集额外支持 GMT。
- 所有模块级表必须包含唯一 `module_uid`：`<network>__<method>__<subtype>__<module>`。
- manifest 中的文件路径使用 `/` 和相对路径，并相对于该 manifest 所在目录解析。
- 每一步既返回 R data.frame，也把主输出路径写入返回对象属性。
- 输入和输出路径均由函数参数或 YAML 指定，不硬编码盘符。
- `strict = TRUE` 时，身份冲突、覆盖不全、非有限值和缺失交接文件立即报错。
- 新输出目录必须不存在，避免静默覆盖正式结果。

## 2. 接口总表

| 步骤 | 函数 | 主输入 | 主输出 |
|---|---|---|---|
| 04 | `module_selection()` | ModuleDivision 根目录 | `module_manifest.tsv` |
| 05 | `functional_annotation()` | `module_manifest.tsv`、基因集 | `module_annotation.tsv` |
| 06 | `drug_response_analysis()` | manifest、药敏根目录或索引 | `drug_response_summary.tsv` |
| 06A | `prepare_module_features()` | manifest、Core34 表 | `module_features.tsv` |
| 06B | `run_nested_ml_scoring()` | features、feature schema | `ml_scores.tsv`、`ml_top10.tsv` |
| 06B | `prepare_ml_scores()` | features、已有概率表 | 同上 |
| 06C | `triage_modules()` | manifest、ML Top10、生存、药敏汇总 | `selected_modules.tsv` |
| 07 | `run_SEQCre()` | selected、蛋白序列、药物 SMILES | `seq_smiles_manifest.tsv` |
| 08 | `predict_BA()` | selected、sequence/SMILES、结合分数/函数 | `binding_scores.tsv` |
| 09 | `process_prs_dti()` | selected、binding、敏感性表/函数 | `final_candidates.tsv` |

第 1–3 步尚未在本包实现，其 ModuleDivision 输入格式在第 3 节说明。

## 3. 第 4 步：ModuleSelection

### 3.1 输入

`base_input_path` 指向 subnetDR 第 3 步 ModuleDivision 根目录。每个亚型/网络/方法需要节点-模块表和带模块标签的边表。标准最小列：

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

模块 0 表示跨模块或未归属边，不进入模块内部边文件。筛选规则为 `module_size > numberCutoff`。

### 3.2 `module_manifest.tsv`

```text
module_uid, legacy_module_id, network, method, subtype, module,
module_size, edge_count, self_loops_excluded,
node_file, edge_file, prefilter_pass, prefilter_reason,
node_sha256, edge_sha256,
source_module_file, source_node_file, source_edge_file
```

每个 `node_file` 为 `node,module`，每个 `edge_file` 为 `node1,node2,module`。边文件只含本模块内部边，并排除自环。

## 4. 第 5 步：ModuleAnnotation

### 4.1 自定义基因集输入

长表最小列：

```text
term_id  gene
GO:1     TP53
GO:1     MDM2
```

可选列：`database,description`。GMT 每行使用 `term_id,description,gene1,gene2,...`。背景基因文件可为单列 gene 表；也可由 `background_genes` 直接传入字符向量。

### 4.2 `module_annotation.tsv`

```text
module_uid, network, method, subtype, module,
database, category, term_id, description,
gene_ratio, background_ratio, rich_factor, fold_enrichment,
p_value, p_adjust, q_value, gene_count, gene_ids,
analysis_status, status_reason
```

具体列以输出表为准。无显著条目的模块仍保留状态行，避免被误认为未运行。

其他输出：

- `module_annotation_top<N>.tsv`
- `module_annotation_qc.tsv`
- `module_annotation_parameters.tsv`
- 可选的逐模块注释文件

## 5. 第 6 步：DrugResponse

### 5.1 输入方式 A：按面板命名的根目录

```r
c(
  PRISM = "path/to/PRISM",
  GDSC1 = "path/to/GDSC1",
  GDSC2 = "path/to/GDSC2",
  CTRP2 = "path/to/CTRP2"
)
```

目录结构可沿用 subnetDR 的网络/方法/亚型/模块布局。解析后立即固化为显式来源索引。

### 5.2 输入方式 B：显式索引

必需列：

```text
module_uid, drug_panel, drn_file, drn_info_file
```

可选列：

```text
drug_level_file, prediction_file
```

两种输入方式必须二选一。

### 5.3 标准 DRN

DRN 至少提供可标准化为以下语义的字段：

```text
protein, drug, p_value
```

可选 `p_adjust,effect_size,effect_direction,tested_samples`。DRN-info 至少需要：

```text
node, type
```

其中 `type` 为 `protein` 或 `drug`。

### 5.4 `drug_response_summary.tsv`

```text
module_uid, network, method, subtype, module, module_size,
drug_panel, drug_number, tested_drug_number, drug_response_density,
drn_drug_number, drn_edge_number, p_adjust_method, significance_cutoff,
drn_file, drn_info_file, analysis_status, status_reason
```

`drug_response_density = drug_number / module_size`。

其他输出：

- `drug_response_hits.tsv`
- `drug_response_source_index.tsv`
- `drug_response_coverage.tsv`
- 标准化的逐模块 DRN/DRN-info

## 6. 第 6A 步：ModuleFeatures

### 6.1 特征输入

每个模块一行，需包含 `module_uid`；或包含能唯一映射到 manifest 的模块身份列。特征列必须为数值且有限。推荐：

```text
module_uid, module_size, <34 feature columns>
```

如果未显式给出 `feature_columns`，程序在排除身份和 QC 列后选择特征，并要求数量等于 `expected_feature_count`。

### 6.2 `module_features.tsv`

```text
module_uid, network, method, subtype, module,
<按 feature_schema.json 固定顺序的 34 个特征>,
feature_missing_count, feature_qc_pass, feature_qc_reason
```

`module_size` 可以作为 Core34 的一个已声明特征；身份列不进入模型。

其他输出：

- `feature_schema.json`：特征版本、数量、顺序、来源路径和 SHA256
- `module_features_qc.tsv`：映射、缺失和模块大小一致性

## 7. 第 6B 步：MLScoring

### 7.1 已有概率表输入

```text
module_uid, prob_C1, prob_C2, prob_C3, prob_C4
```

允许附带身份、预测亚型和模型元数据，但四个概率和全部 manifest 模块必须完整。每行概率和应为 1（允许数值误差）。

### 7.2 内置嵌套 OOF 输入

- 第 6A 步的 `module_features.tsv`
- 同目录的 `feature_schema.json`
- outer/inner CV、随机种子、Top-K 和 `n_jobs`

`fast` 使用论文已选 Gradient Boosting 参数；`paper` 在每个外层训练折中执行完整内层网格。

### 7.3 `ml_scores.tsv`

```text
module_uid, network, method, subtype, module,
true_subtype, predicted_subtype,
prob_C1, prob_C2, prob_C3, prob_C4,
target_subtype_probability, probability_margin,
rank_in_subtype, score_type, model_version
```

### 7.4 `ml_top10.tsv`

包含 `ml_scores.tsv` 全部列，另加：

```text
target_subtype, top_k, ml_gate, ml_gate_reason
```

文件名实际随 `top_k` 变化，例如 `top_k=10` 时为 `ml_top10.tsv`。

其他输出：

- `ml_gate_audit.tsv` 或 Python fold/metadata 文件
- `fitted_model.joblib`（内置 Python 路径）

## 8. 第 6C 步：ModuleTriage

### 8.1 生存输入

每个 ML Top-K 模块一行。标准列：

```text
module_uid, n_samples, n_events, matched_gene_fraction,
cox_hr_per_sd, cox_p, cox_fdr, cutpoint,
logrank_p, logrank_fdr, hr_high_vs_low, survival_direction
```

常用同义列会被标准化。`survival_direction` 默认要求 `High_score_worse`。

### 8.2 `module_evidence.tsv`

把 ML Top-K、manifest、生存和主药敏面板按 `module_uid` 一对一连接，并添加：

```text
prognosis_gate, prognosis_gate_reason,
module_size_gate, module_size_gate_reason,
drug_response_gate, drug_response_gate_reason, all_gates_pass
```

### 8.3 `module_filtering_stepwise.tsv`

```text
module_uid, subtype, gate_order, gate_name,
gate_pass, cumulative_pass, gate_reason
```

### 8.4 `selected_modules.tsv`

```text
module_uid, network, method, subtype, module, module_size,
primary_drug_panel, node_file, edge_file, drn_file, drn_info_file,
target_subtype_probability, probability_margin, ml_rank_in_subtype,
survival_direction, logrank_p, logrank_fdr, hr_high_vs_low,
drug_number, tested_drug_number, drug_response_density, drn_edge_number,
selection_rank, selection_reason
```

默认每个亚型恰好一行。四个文件路径均指向第 6C 输出目录内复制的自包含交接文件。

## 9. 第 7 步：SEQCre

### 9.1 查找表输入

蛋白序列：

```text
node, sequence
```

药物结构：

```text
node, smiles
```

可识别 protein/target/gene 和 drug/compound 等常见 ID 同义列。

### 9.2 `seq_smiles_manifest.tsv`

```text
module_uid, subtype, sequence_file, smiles_file,
required_protein_number, matched_protein_number,
required_drug_number, matched_drug_number,
analysis_status, status_reason
```

每模块的 `protein_sequences.tsv` 为 `node,sequence,sequence_status`；`drug_smiles.tsv` 为 `node,SMILES,smiles_status`。

## 10. 第 8 步：BindingScore

### 10.1 结合分数输入

```text
module_uid, target_name, drug_name, binding_score
```

每个 selected DRN 药物-靶点对必须唯一。也可传入 `predictor(dpi)`，其中 DPI 包含：

```text
module_uid, target_name, drug_name, target_seq, drug_smiles
```

### 10.2 `binding_scores.tsv`

```text
module_uid, target_name, drug_name, target_seq, drug_smiles,
binding_score, binding_rank, score_direction
```

另输出 `binding_score_manifest.tsv` 和逐模块 DPI/结合分数文件。

## 11. 第 9 步：PScore

### 11.1 靶点敏感性输入

```text
module_uid, target_name, sensitivity
```

也可传入 `sensitivity_function(edges, selected_module_row)`。函数返回 target 和 sensitivity 列。

### 11.2 `perturbation_scores.tsv`

包含 selected 模块身份、步骤 8 的全部结合字段，以及：

```text
sensitivity, perturbation_score, perturbation_rank_in_module
```

其中：

```text
perturbation_score = binding_score * sensitivity
```

### 11.3 `final_candidates.tsv`

是每亚型按扰动分数排序后保留的前 `top_n` 行，在扰动分数表基础上增加：

```text
final_rank_in_subtype
```

另输出标准化 `target_sensitivity.tsv`。

## 12. 总运行器配置要求

`run_ML_SnpDR()` 从配置中的 `paths` 读取外部输入和各阶段输出目录。关键外部输入字段为：

```text
module_division_dir
annotation_gene_set_file
annotation_background_gene_file
drug_response_input_index_file 或 drug_response.input_roots
module_feature_input_file
ml_score_input_file（为空时运行内置嵌套 OOF）
survival_input_file
protein_sequence_file
drug_smiles_file
binding_score_input_file
sensitivity_input_file
```

当 `from` 指向中间步骤时，配置必须指向已存在的上游主输出。`dry_run=TRUE` 只检查配置和阶段区间，不创建科学结果。
