# 中文使用说明

按 Step 1–5 完成 BAM 处理、特征矩阵构建、集成模型训练、验证和盲测。训练默认使用 `marker_file/` 中的研究 marker，癌症检测和组织溯源分别读取各自的 7 份列表。

## 运行模拟示例

按照[安装说明](ENVIRONMENT.md)配置环境，在仓库根目录运行：

```bash
Rscript scripts/check_environment.R all
python3 examples/run_demo.py --out work/demo --threads 1
```

示例包含一个模拟 BAM，以及 Normal、肺癌、结直肠癌、肝癌四组模拟数据：训练集 160 例、验证集 40 例、盲测集 40 例。程序会生成覆盖完整 marker 的模拟矩阵，训练全部五种算法和两层模型，并预测模拟 BAM。完成后生成 `work/demo/COMPLETED.txt`。每次运行使用新的输出目录。

原始研究模型和参与者数据不包含在仓库中。示例运行时会重新拟合模型；示例指标仅用于检查程序。

## 五个步骤

| 步骤 | 脚本 | 输入与输出 |
| --- | --- | --- |
| Step 1 | `scripts/step1_bam_to_features.py` | BAM、hg19 CpG 注释和 BED → mHap、CGI 统计、EM/BPM |
| Step 2 | `scripts/step2_build_matrices.R` | 样本 profile、七类特征表，或 CGI 加 EM/BPM 矩阵 → 七类 R 矩阵 |
| Step 3 | `scripts/step3_train_ensemble.R` | 训练矩阵和默认 marker → 新模型、折分、折外预测、阈值 |
| Step 4 | `scripts/step4_predict.R` | 新模型和验证/盲测矩阵 → 检测分数、阳性样本溯源概率 |
| Step 5 | `scripts/step5_evaluate.R` | 固定预测和独立标签表 → 检测与溯源统计 |

每一步的命令见[README](../README.md)。甲基化特征为 MM、PDR、CHALM、MHL、Entropy，片段组学特征为 EM 和 BPM。第一层使用 LR、GBM、线性 SVM、XGBoost、CatBoost，第二层使用线性 SVM。两个任务均为五折，共生成 360 个新的模型文件。

`--preset demo` 缩短调参和 boosting 迭代；`--preset study` 使用原脚本设置。两者都使用完整 marker、五种算法和两层结构。

## 已有 mHap 和 EM/BPM 矩阵

填写 `examples/processed_cohorts_template.tsv`。每行一个参与者样本，指定训练、验证或盲测分组；每行填写 `mhap` 或 `cgi` 之一。已有 CGI 统计时可以跳过 mHap 计算。

```bash
python3 workflow/processed_cohorts.py \
  --manifest /path/to/cohorts.tsv \
  --em /path/to/EM.RData --bpm /path/to/BPM.RData \
  --motif-format counts --cpg /path/to/hg19_CpG.gz \
  --out work/cohorts --jobs 2

Rscript scripts/step3_train_ensemble.R \
  --train work/cohorts/matrices/train --out work/model \
  --preset study --threads 4
```

EM/BPM 支持 TSV、CSV、RDS 和单对象 RData，可以是样本×motif，也可以是 motif×样本。所有 4,096 个标准六碱基 motif 都需提供。默认输入是原始计数；频率输入需用 `--motif-format` 指明归一化方式。

只有标准 motif 的计数或频率已经丢弃了含 N 的信息时，需在样本表提供 `EM_valid_fraction` 和 `BPM_valid_fraction`，即标准 motif 计数占原总计数的比例。盲测会据此使用原分析的分母。具体格式见[输入说明](INPUTS.md)。

## marker

`marker_file/` 保存 14 份 CSV，`markers.tsv` 为合并后的训练输入。五类甲基化 marker 使用 hg19 的 `chr.start.end` 坐标；EM/BPM 使用六碱基序列。计数见[marker 说明](../marker_file/README.md)。

Step 3 默认读取该表，训练前检查所有 marker 是否存在于矩阵中，并将实际使用的列表保存为 `markers_used.tsv`。验证和盲测使用训练时保存的特征顺序。小型模拟种子表中的 marker 仅用于生成初始示例，不是默认训练列表。

验证和盲测的诊断标签放在单独的 `truth.tsv` 中，Step 5 才读取。运行输出保存在 `work/`。
