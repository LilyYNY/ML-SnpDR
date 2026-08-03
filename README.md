# ML-SnpDR

Machine-Learning-Guided Subtype-Specific Network Drug Repurposing

ML-SnpDR 基于 [subnetDR](https://github.com/LilyYNY/subnetDR) 重构 LUAD 亚型特异性网络药物重定位流程。核心改变发生在原流程第 4–6 步之后：先对第 4 步保留的全部模块完成注释、药物响应、Core34 特征和机器学习评分；每个亚型取机器学习 Top10，再结合生存和药物响应证据选出一个最优模块；最后只让这些入选模块进入 subnetDR 的第 7–9 步。

当前第 1–9 步均已实现。默认总入口从表达矩阵和亚型表开始，依次生成差异表达、亚型 PPI、模块划分、模块清单、注释、药物响应、Core34、ML Top10、最优模块和最终药物-靶点候选。也可以从任意已存在的中间结果开始。

## 流程

```mermaid
flowchart TD
    E["expression + subtype phenotype"] --> S1["01 DEPs：Wilcoxon + BH"]
    S1 --> D["differential_expression.tsv"]
    D --> S2["02 NetworkConstruction：亚型上调蛋白诱导 PPI"]
    S2 --> N["network_manifest.tsv"]
    N --> S3["03 ModuleDivision：Louvain / WF"]
    S3 --> S4["04 ModuleSelection"]
    S4 --> M["module_manifest.tsv：全部大小预筛模块"]
    M --> S5["05 ModuleAnnotation：全部模块"]
    M --> S6["06 DrugResponse：全部模块 × 全部面板"]
    M --> S6A["06A ModuleFeatures：全部模块 Core34"]
    S6A --> S6B["06B MLScoring：全部模块概率 + 每亚型 Top10"]
    S6B --> S6C["06C ModuleTriage：Top10 生存 + 药物证据筛选"]
    S6 --> S6C
    M --> S6C
    S6C --> SEL["selected_modules.tsv：每亚型最优模块"]
    SEL --> S7["07 SEQCre：仅入选模块"]
    S7 --> S8["08 BindingScore：仅入选模块"]
    S8 --> S9["09 PScore：仅入选模块"]
```

机器学习概率表示模块的亚型代表性，不等同于疗效或可成药性。因此，ML Top10 必须继续经过生存、模块大小和药物响应证据筛选。

## 功能状态

| 步骤 | 函数 | 分析范围 | 主输出 | 状态 |
|---|---|---|---|---|
| 01 | `run_diff_expr_analysis()` | 全部样本；亚型 vs 其余样本 | `differential_expression.tsv` | 可用 |
| 02 | `run_network_construction()` | 每亚型显著上调蛋白 × PPI 来源 | `network_manifest.tsv`、亚型 PPI | 可用 |
| 03 | `module_division()` / `subtype_module()` | 亚型 × 网络 × Louvain/WF | `module_division_manifest.tsv` | 可用 |
| 04 | `module_selection()` | 全部大小预筛模块 | `module_manifest.tsv` | 可用 |
| 05 | `functional_annotation()` | manifest 全部模块 | `module_annotation.tsv` | 可用 |
| 06 | `drug_response_analysis()` | manifest 全部模块 × 面板 | `drug_response_summary.tsv`、DRN | 可用 |
| 06A | `prepare_module_features()` | manifest 全部模块 | `module_features.tsv`、`feature_schema.json` | 可用 |
| 06B | `prepare_ml_scores()` / `run_nested_ml_scoring()` | 全部模块；每亚型 Top10 | `ml_scores.tsv`、`ml_top10.tsv` | 可用 |
| 06C | `triage_modules()` | 每亚型 ML Top10 | `selected_modules.tsv`、筛选轨迹 | 可用 |
| 07 | `run_SEQCre()` | 每亚型最优模块 | `seq_smiles_manifest.tsv` | 可用 |
| 08 | `predict_BA()` | 每亚型最优模块 | `binding_scores.tsv` | 可用 |
| 09 | `process_prs_dti()` | 每亚型最优模块 | `perturbation_scores.tsv`、`final_candidates.tsv` | 可用 |

此外还提供 YAML 配置合并与校验、网络/算法名称标准化、`module_uid` 创建和解析、阶段注册表、dry-run、输入覆盖度检查、SHA256 来源追踪和逐步 QC。

## 每一步输入输出总览

下表用于快速准备文件。完整字段、允许的同义列、状态列和严格校验规则见 [每一步输入输出格式](docs/io-contracts.md)。

| 步骤 | 自定义输入 | 最小输入格式 | 主要输出 | 下游用途 |
|---|---|---|---|---|
| 01 DEPs | 表达矩阵、样本亚型表、检出率/FC/P 阈值 | 表达宽表：`gene,<sample...>`；表型：`Sample,Subtype` | `differential_expression.tsv`、显著结果、统计汇总、兼容 Excel | 第 2 步亚型蛋白集合 |
| 02 NetworkConstruction | 差异表达表、PPI 索引或命名边文件、分数阈值 | PPI 索引至少：`network,edge_file`；边表至少两端点列 | `network_manifest.tsv`、每亚型/网络 `ppi_<subtype>.txt` | 第 3 步网络输入 |
| 03 ModuleDivision | network manifest、网络/算法、随机种子、WF 阈值 | manifest：`network,subtype,ppi_file` | `module_division_manifest.tsv`、节点/边模块文件 | 第 4 步模块筛选输入 |
| 04 ModuleSelection | ModuleDivision 根目录、亚型、网络、算法、大小阈值 | 节点表：`node,module`；边表：`node1,node2,module` | `module_manifest.tsv`、逐模块 `nodes.tsv/edges.tsv` | 第 5、6、6A、6C 步共享的模块全集 |
| 05 ModuleAnnotation | `module_manifest.tsv`、GMT/TSV/CSV/XLSX 基因集、背景基因 | 长表：`term_id,gene`，可选 `database,description` | `module_annotation.tsv`、Top-N、QC | 模块生物学解释和报告 |
| 06 DrugResponse | manifest；四个面板根目录或显式索引 | 索引：`module_uid,drug_panel,drn_file,drn_info_file` | `drug_response_summary.tsv`、`drug_response_hits.tsv`、标准化 DRN | 第 6C 步药物响应证据 |
| 06A ModuleFeatures | manifest、全模块 Core34 表 | `module_uid` + 34 个数值特征 | `module_features.tsv`、`feature_schema.json`、QC | 第 6B 步机器学习输入 |
| 06B MLScoring | features + schema，或已有全模块概率表 | 概率表：`module_uid,prob_C1,prob_C2,prob_C3,prob_C4` | `ml_scores.tsv`、`ml_top10.tsv` | 第 6C 步每亚型候选集合 |
| 06C ModuleTriage | manifest、ML Top10、生存表、药物响应汇总 | 生存表至少需要 `module_uid,logrank_p,survival_direction` | `module_evidence.tsv`、`module_filtering_stepwise.tsv`、`selected_modules.tsv` | 建立第 7–9 步的入选模块边界 |
| 07 SEQCre | selected、蛋白序列表、药物 SMILES 表 | `node,sequence`；`node,smiles` | `seq_smiles_manifest.tsv`、逐模块序列/SMILES | 第 8 步 DPI 和结合预测输入 |
| 08 BindingScore | selected、sequence/SMILES manifest、结合分数表或预测函数 | `module_uid,target_name,drug_name,binding_score` | `binding_scores.tsv`、`binding_score_manifest.tsv` | 第 9 步结合证据 |
| 09 PScore | selected、binding、敏感性表或 ENM/PRS 函数 | `module_uid,target_name,sensitivity` | `perturbation_scores.tsv`、`final_candidates.tsv` | 最终每亚型药物-靶点候选 |

这里有两个不能混淆的数据边界：

- `module_manifest.tsv`：第 4 步大小预筛后的全部模块；第 5、6、6A 和 6B 步不得提前删除其中的模块。
- `selected_modules.tsv`：第 6C 步每亚型筛出的最优模块；第 7、8、9 步只能读取这份表列出的模块和文件。

## 安装

```bash
git clone https://github.com/LilyYNY/ML-SnpDR.git
cd ML-SnpDR
conda env create -f environment.yml
conda activate mlsnpdr
R CMD INSTALL .
```

也可以分别安装：

```r
remotes::install_github("LilyYNY/ML-SnpDR")
```

```bash
python -m pip install -e ".[dev]"
```

## 统一模块主键

全部模块级表使用：

```text
<network_slug>__<method_slug>__<subtype>__<module>
```

例如：

```r
make_module_uid("physicalPPIN", "Louvain", "C3", "M10")
# physicalppin__louvain__C3__M10
```

任何一步都不能靠目录名重新猜测模块身份，必须通过 `module_uid` 和显式文件路径传递。

## 单步运行

以下示例展示每一步怎样直接使用上一步的输出。每个函数都允许自定义输入和输出位置，且不调用全局 `setwd()`。

### 第 1 步：亚型 vs 其余样本差异表达

**方法**：与 subnetDR 一致，对每个亚型分别进行 one-versus-rest Wilcoxon 秩和检验，使用 BH 校正；蛋白必须在目标组和其余样本组同时满足检出率阈值，再根据校正 P 值和 fold change 标记为 `up`、`down` 或 `non_significant`。

**输入**：表达矩阵、样本亚型表，以及可自定义的检出率、fold change、P 值和多重校正参数。表达输入是宽表，第一列是基因/蛋白，其余列是样本：

```text
gene    S001    S002    S003
EGFR    10.2    8.7     4.1
TP53    3.2     3.8     5.6
```

亚型输入至少包含：

```text
Sample  Subtype
S001    C1
S002    C1
S003    C2
```

```r
deps <- run_diff_expr_analysis(
  expression_file = "data/raw/expression.xlsx",
  phenotype_file = "data/raw/subtype.xlsx",
  output_dir = "results/01_differential_expression",
  gene_column = "gene",
  sample_column = "Sample",
  subtype_column = "Subtype",
  detection_threshold = 0.80,
  fc_threshold = 2,
  p_threshold = 0.01,
  p_adjust_method = "BH",
  pseudocount = 1,
  write_legacy_excel = TRUE
)

diff_file <- attr(deps, "result_file")
```

**输出**：

- `differential_expression.tsv`：所有亚型、所有蛋白的长表，是第 2 步主输入；
- `differential_expression_significant.tsv`：仅 `up/down`；
- `differential_expression_summary.tsv`：每亚型各标签数量；
- `diff_expression_results_all.xlsx`：兼容 subnetDR 的逐亚型工作表。

标准长表包含 `subtype,gene,mean_target,mean_other,fold_change,log2_fold_change,p_value,p_adjust,detection_target,detection_other,label,n_target,n_other`。

### 第 2 步：上调蛋白构建亚型特异 PPI

**方法**：默认使用第 1 步每个亚型标记为 `up` 的蛋白，并在 String、physicalPPIN、chengF 中保留两个端点都属于该集合的相互作用。PPI 来源不再写死在 `./PPI`，而由显式索引指定。

**输入**：第 1 步的 `differential_expression.tsv`，以及 PPI 索引或命名 PPI 边文件。最简单的 PPI 索引为：

```text
network        edge_file
String         PPI/String/string_symbol_edges.tsv
physicalPPIN   PPI/PhysicalPPIN/physicalppi.tsv
chengF         PPI/ChengF/chengf_symbol_edges.tsv
```

若原始 PPI 使用数据库 ID，可增加以下列：

```text
network, edge_file, node1_column, node2_column,
score_column, score_min, delimiter,
mapping_file, mapping_id_column, mapping_symbol_column, mapping_delimiter
```

`edge_file` 和 `mapping_file` 的相对路径均相对于 PPI 索引文件解析。String 空格分隔文件可写 `delimiter=whitespace`；映射表为制表符时写 `mapping_delimiter=tab`。

```r
networks <- run_network_construction(
  diff_file = diff_file,
  ppi_index_file = "data/raw/ppi_index.tsv",
  output_dir = "results/02_network_construction",
  ppi_method = c("String", "physicalPPIN", "chengF"),
  ppiScore = 400,
  include_labels = "up",
  subtypes = c("C1", "C2", "C3", "C4")
)

network_manifest_file <- attr(networks, "manifest_file")
```

**输出**：

- `network_manifest.tsv`：一行一个亚型-网络组合；
- `<subtype>/<network>/ppi_<subtype>.txt`：标准 `node1,node2` 边表；
- 网络节点数、边数、来源路径、状态和 SHA256。

### 第 3 步：Louvain/WF 模块划分

**方法**：`Louvain` 直接在无向、去重、无自环 PPI 上运行。`WF` 分别运行 edge-betweenness 和 label propagation，对两组模块做超几何重叠检验，再按 P 值和重叠大小确定互不重叠的共识模块。

**输入**：第 2 步的 `network_manifest.tsv`，以及需要运行的亚型、网络、模块算法、随机种子和 WF 显著性阈值。

```r
divisions <- module_division(
  network_manifest_file = network_manifest_file,
  output_base_path = "results/03_module_division",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF"),
  subtypes = c("C1", "C2", "C3", "C4"),
  seed = 123,
  wf_pvalue_cutoff = 0.05,
  strict = TRUE
)

module_division_dir <- attr(divisions, "output_dir")
```

**输出**：每个亚型-网络-算法组合生成：

- `node_Module_<network>_<method>.txt`：`node,module`；
- `edges_<network>_<method>.txt`：`node1,node2,module`，跨模块边为 0；
- `edge_Module_<network>_<method>.txt`：仅模块内部边；
- 根目录 `module_division_manifest.tsv`：记录全部组合、计数、路径、参数和 SHA256。

### 第 4 步：ModuleSelection → manifest

**方法**：统计第 3 步每个模块的唯一节点数，保留 `module_size > numberCutoff` 的全部模块；排除模块 0 和自环，同时保留 subnetDR 聚合文件，并为每个模块生成唯一 `module_uid`、独立节点/边文件和 SHA256。

**输入**：第 3 步的模块划分根目录、亚型表、网络名、模块算法名和模块大小阈值；目录内节点表至少包含 `node,module`，边表至少包含 `node1,node2,module`。

```r
library(MLSnpDR)

modules <- module_selection(
  subtype_file = "data/subtype.xlsx",
  base_input_path = module_division_dir,
  base_output_path = "results/04_module_selection",
  network_method = c("String", "physicalPPIN", "chengF"),
  module_method = c("Louvain", "WF"),
  numberCutoff = 9,
  strict = TRUE
)

manifest_file <- attr(modules, "manifest_file")
```

**输出**：

- `standardized/module_manifest.tsv`
- `standardized/modules/<module_uid>/nodes.tsv`
- `standardized/modules/<module_uid>/edges.tsv`
- subnetDR 兼容的聚合节点/边文件和 QC

如果已经有 ModuleSelection 输出，可使用 `build_module_manifest()` 直接适配，无需重跑模块筛选。

### 第 5 步：全部模块功能注释

**方法**：以模块节点为查询集、背景基因为 universe，对每个基因集执行超几何过度富集检验，并按指定方法做多重检验校正。注释只用于解释，不会改变进入第 6A/6B 步的模块全集。

**输入**：第 4 步的 `module_manifest.tsv`、GMT/TSV/CSV/XLSX 基因集、背景基因及富集参数。长表基因集至少包含 `term_id,gene`，可选 `database,description`。

```r
annotation <- functional_annotation(
  module_manifest_file = manifest_file,
  base_output_path = "results/05_module_annotation",
  gene_set_file = "data/gene_sets.gmt",
  background_gene_file = "data/background_genes.tsv",
  pAdjustMethod = "BH",
  pAdjustCutoff = 0.05,
  top_n = 15
)
```

`gene_set_file` 可为 TSV、CSV、XLSX 或 GMT。若不提供自定义基因集，可通过 `databases` 调用 `msigdbr`。

**输出**：`module_annotation.tsv`、`module_annotation_top<N>.tsv`、`module_annotation_qc.tsv`、参数表和可选逐模块注释文件。

### 第 6 步：全部模块药物响应标准化

**方法**：对 `module_manifest.tsv` 的全部模块和全部指定药敏面板建立明确的模块-面板索引，统一已有 DRN、DRN-info、药物显著性和效应方向，并执行全覆盖审计。该步骤不根据药敏结果提前删除模块。

**输入**：第 4 步的 `module_manifest.tsv`，以及 PRISM、GDSC1、GDSC2、CTRP2 等药物响应根目录或显式索引；索引至少包含 `module_uid,drug_panel,drn_file,drn_info_file`。

```r
drug <- drug_response_analysis(
  module_manifest_file = manifest_file,
  drug_response_path = "results/06_drug_response",
  drug_response_roots = c(
    PRISM = "existing/PRISM",
    GDSC1 = "existing/GDSC1",
    GDSC2 = "existing/GDSC2",
    CTRP2 = "existing/CTRP2"
  ),
  panels = c("PRISM", "GDSC1", "GDSC2", "CTRP2"),
  strict = TRUE
)

drug_summary_file <- attr(drug, "summary_file")
```

也可用 `drug_response_index_file` 代替多个根目录。索引至少包含：

```text
module_uid, drug_panel, drn_file, drn_info_file
```

可选列为 `drug_level_file`、`prediction_file`。第 6 步必须覆盖 manifest 中的全部模块，而不是先按药物结果删模块。

**输出**：`drug_response_summary.tsv`、`drug_response_hits.tsv`、`drug_response_source_index.tsv`、`drug_response_coverage.tsv` 和标准化逐模块 DRN/DRN-info。

### 第 6A 步：全部模块 Core34 特征

**方法**：将 Core34 特征表与第 4 步 manifest 做一对一身份映射，固定特征顺序，检查全模块覆盖、数值有限性、缺失数和模块大小一致性，并记录来源哈希。

**输入**：第 4 步的 `module_manifest.tsv` 和全模块 Core34 表。特征表必须一行一个模块，包含 `module_uid`，或包含能够唯一构造身份的 network/method/subtype/module 列，以及 34 个数值特征。

```r
features <- prepare_module_features(
  module_manifest_file = manifest_file,
  feature_file = "data/all_module_core34.tsv",
  output_dir = "results/06A_module_features",
  expected_feature_count = 34,
  feature_set = "core34_v1",
  strict = TRUE
)

feature_file <- attr(features, "feature_file")
feature_schema_file <- attr(features, "schema_file")
```

输出固定特征顺序，并检查覆盖度、缺失值、有限值和模块大小一致性。

**输出**：`module_features.tsv`、`feature_schema.json` 和 `module_features_qc.tsv`。

### 第 6B 步：嵌套 OOF 评分和每亚型 Top10

**方法**：在全部模块上运行重复嵌套分层交叉验证 Gradient Boosting，使用外层 OOF 概率避免训练内评分，并按模块真实亚型的目标概率生成亚型内排名；通过预测亚型匹配等 gate 后每亚型保留 Top10。也可以导入满足同一概率契约的外部模型结果。

**输入**：第 6A 步的 `module_features.tsv` 和 `feature_schema.json`；若导入已有模型结果，则概率表必须包含 `module_uid,prob_C1,prob_C2,prob_C3,prob_C4`。

直接运行内置 Python 嵌套 OOF：

```r
scores <- run_nested_ml_scoring(
  module_features_file = feature_file,
  feature_schema_file = feature_schema_file,
  output_dir = "results/06B_ml_scoring",
  mode = "paper",
  outer_splits = 5,
  outer_repeats = 5,
  inner_splits = 3,
  top_k = 10
)

ml_top_file <- attr(scores, "top_file")
```

也可以导入已有概率表：

```r
scores <- prepare_ml_scores(
  module_features_file = feature_file,
  score_file = "data/all_module_probabilities.tsv",
  output_dir = "results/06B_ml_scoring",
  top_k = 10,
  require_predicted_subtype_match = TRUE
)
```

程序校验全模块覆盖、概率和、预测亚型、目标亚型概率和亚型内排名。

**输出**：`ml_scores.tsv`、`ml_top10.tsv`、gate 审计、CV fold/模型元数据和 `fitted_model.joblib`（内置 Python 路径）。

### 第 6C 步：Top10 + 生存 + 药物响应 → 每亚型最优模块

**方法**：只在每亚型 ML Top10 内连接模块生存结果和第 6 步主药敏面板证据，依次执行预后方向、log-rank、模块大小和药物响应 gate，再按预设证据优先级排序。

**输入**：第 4 步的 `module_manifest.tsv`、第 6B 步的 `ml_top10.tsv`、模块生存表和第 6 步的 `drug_response_summary.tsv`。生存表至少包含 `module_uid,logrank_p,survival_direction`。

```r
selected <- triage_modules(
  module_manifest_file = manifest_file,
  ml_top_file = ml_top_file,
  survival_file = "data/module_survival.tsv",
  drug_response_summary_file = drug_summary_file,
  output_dir = "results/06C_module_triage",
  primary_drug_panel = "PRISM",
  min_module_size = 30,
  required_direction = "High_score_worse",
  logrank_p_max = 0.05,
  select_n_per_subtype = 1
)

selected_file <- attr(selected, "selected_file")
```

默认筛选顺序为：ML 合格 Top10 → `High_score_worse` 且 log-rank P ≤ 0.05 → 模块大小 ≥ 30 → 有主药敏面板证据 → 按显著药物数、药物响应密度和目标亚型概率排序 → 每个亚型选择 1 个模块。

输出目录会复制入选模块的 node、edge、DRN 和 DRN-info 文件，因此第 7 步只依赖 `selected_modules.tsv` 及其相对路径，不再扫描全部模块目录。

**输出**：`module_survival.tsv`、`module_evidence.tsv`、`module_filtering_stepwise.tsv`、`selected_modules.tsv` 及入选模块的四类自包含交接文件。

### 第 7 步：仅为入选模块生成 sequence/SMILES

**方法**：逐行读取 `selected_modules.tsv` 的 DRN-info，只提取实际出现的 protein/drug，然后与显式序列和 SMILES 查找表连接并检查覆盖度，不递归扫描其他模块。

**输入**：第 6C 步的 `selected_modules.tsv`、蛋白序列表 `node,sequence` 和药物 SMILES 表 `node,smiles`。

```r
seqs <- run_SEQCre(
  input_base = selected_file,
  output_base = "results/07_sequence_smiles",
  protein_sequence_file = "data/protein_sequences.tsv",
  drug_smiles_file = "data/drug_smiles.tsv"
)

seq_manifest_file <- attr(seqs, "manifest_file")
```

常见同义列会被标准化，但推荐使用上述四个列名。

**输出**：`seq_smiles_manifest.tsv`，以及每个入选模块的 `protein_sequences.tsv` 和 `drug_smiles.tsv`。

### 第 8 步：仅为入选模块计算结合分数

**方法**：从入选 DRN 构造 drug-target pair，连接第 7 步的 target sequence 和 drug SMILES，再导入已有结合分数或调用用户预测函数；根据 `lower_better/higher_better` 生成模块内排名。

**输入**：第 6C 步的 `selected_modules.tsv`、第 7 步的 `seq_smiles_manifest.tsv`，以及结合分数表或用户预测函数。结合分数表包含 `module_uid,target_name,drug_name,binding_score`。

```r
binding <- predict_BA(
  selected_modules_file = selected_file,
  seq_smiles_manifest_file = seq_manifest_file,
  output_base = "results/08_binding_score",
  binding_score_file = "data/binding_scores.tsv",
  score_direction = "lower_better"
)

binding_file <- attr(binding, "scores_file")
```

也可通过 `predictor` 传入用户自己的预测函数；文件和函数必须二选一。

**输出**：`binding_scores.tsv`、`binding_score_manifest.tsv` 和逐模块 DPI/结合分数文件。

### 第 9 步：PRS/扰动评分和最终候选

**方法**：从入选模块边文件导入或计算 target sensitivity，与第 8 步结合分数一对一连接，计算 perturbation score，并分别生成模块内和亚型内排名。

**输入**：第 6C 步的 `selected_modules.tsv`、第 8 步的 `binding_scores.tsv`，以及敏感性表或 ENM/PRS 函数。敏感性表包含 `module_uid,target_name,sensitivity`。

```r
prs <- process_prs_dti(
  selected_modules_file = selected_file,
  binding_scores_file = binding_file,
  output_base = "results/09_perturbation_score",
  sensitivity_file = "data/target_sensitivity.tsv",
  top_n = 10
)
```

也可传入 `sensitivity_function`。扰动分数沿用 subnetDR 的定义：

```text
perturbation_score = binding_score × sensitivity
```

**输出**：`target_sensitivity.tsv`、`perturbation_scores.tsv` 和每亚型前 `top_n` 的 `final_candidates.tsv`。

## 一次贯穿第 1–9 步

先在 `config/paper_luad.yml` 中设置原始输入和各阶段输出路径。完整运行至少需要配置：

```text
paths.expression_file
paths.subtype_file
paths.ppi_index_file
paths.module_feature_input_file
paths.survival_input_file
drug_response.input_roots 或 paths.drug_response_input_index_file
paths.protein_sequence_file
paths.drug_smiles_file
paths.binding_score_input_file
paths.sensitivity_input_file
```

`paths.ml_score_input_file` 可选：为空时运行内置嵌套 OOF；提供时导入已有全模块概率。配置完成后只需要一个入口：

```r
library(MLSnpDR)

plan <- run_ML_SnpDR(
  config = "config/paper_luad.yml",
  from = "deps",
  to = "perturbation_score",
  dry_run = TRUE
)

result <- run_ML_SnpDR(
  config = "config/paper_luad.yml",
  from = "deps",
  to = "perturbation_score",
  dry_run = FALSE
)
```

`from="deps"` 是默认值，因此正式运行也可简写为：

```r
result <- run_ML_SnpDR(
  config = "config/paper_luad.yml",
  dry_run = FALSE
)
```

总运行器按阶段注册表顺序执行，并把实际返回的主输出路径直接传给依赖它的下一步。也可只运行连续子区间，例如：

```r
run_ML_SnpDR(
  config = "config/paper_luad.yml",
  from = "module_features",
  to = "module_triage",
  dry_run = FALSE
)
```

若从中间步骤开始，配置中必须提供前置步骤产生的输入文件。配置字段和所有列定义见 [输入输出契约](docs/io-contracts.md)，设计边界见 [流程架构](docs/architecture.md)。

## 当前 LUAD 数据验证

代码已在当前项目数据上完成只读或独立验证目录测试：

- 第 1–3 步：合成表达/PPI 数据完成 DEPs → network manifest → Louvain/WF → ModuleSelection 的连续测试；
- 第 4 步：208 个模块，C1/C2/C3/C4 分别为 63/44/49/52；14,261 条节点记录、95,889 条边，排除 854 个自环。
- 第 6 步：208 模块 × 4 面板 = 832 个模块-面板组合，PRISM/GDSC1/GDSC2/CTRP2 均完整覆盖。
- 第 6A 步：208 个模块 × 34 个特征，缺失模块 0，QC 通过 208。
- 第 6B 步：208 个模块全部评分，每亚型 10 个，共 40 个 ML Top10 模块。
- 第 6C 步：每亚型选择 1 个模块，并确认 4 类交接文件均存在。

这部分结果用于验证代码连接和数据契约，不替代最终论文运行时的固定环境、参数记录和统计复核。

## 测试

```bash
Rscript -e "devtools::test()"
python -m pytest
R CMD check .
```

## 与 subnetDR 的关系

本仓库沿用 subnetDR 的科学工作流和第 7–9 步命名，但移除了这些步骤对全部目录的递归扫描，并把输入范围锁定到 `selected_modules.tsv`。建议同时引用原始 subnetDR 仓库及相应论文/方法来源。

## License

MIT。上游 subnetDR 代码和相关方法的署名与版权归其原作者所有。
