# 中文使用说明

这份代码按 **Step 1–5** 组织：从 BAM 提取 mHap、EM、BPM 和五种甲基化特征，使用指定 marker 训练两层集成模型，再进行验证和盲测。

公开示例包含模拟 BAM，以及肺癌、结直肠癌、肝癌和 Normal 四组模拟数据。训练集 160 例，验证集 40 例，盲测集 40 例，每组人数相同。所有示例编号、性别和分期标签都是模拟的。原始研究的训练模型和参与者数据不在这个仓库中；示例训练时会在本地产生新的模型。

## 一条命令运行完整流程

先按照[依赖安装说明](ENVIRONMENT.md)配置环境，在仓库根目录运行：

```bash
Rscript scripts/check_environment.R all
python3 examples/run_demo.py --out work/demo --threads 1
```

程序会执行 BAM 提取、构建矩阵、训练、预测和评估，并用新训练的模型预测这个模拟 BAM。完成后会生成 `work/demo/COMPLETED.txt`。每次运行应使用一个新的输出目录。

| 步骤 | 脚本 | 输入 → 输出 |
| --- | --- | --- |
| Step 1 | `scripts/step1_bam_to_features.py` | BAM + hg19 CpG 注释 + BED → mHap、甲基化区域统计、EM/BPM |
| Step 2 | `scripts/step2_build_matrices.R` | 每个样本的 profile，或七个 TSV 特征表 → 七类样本×特征矩阵 |
| Step 3 | `scripts/step3_train_ensemble.R` | 训练矩阵 + marker 表 → 新模型、折分、OOF 预测和阈值 |
| Step 4 | `scripts/step4_predict.R` | 新模型 + 未参与训练的矩阵 → 检测结果和阳性样本的溯源概率 |
| Step 5 | `scripts/step5_evaluate.R` | 已固定的预测 + 独立 truth 表 → 检测和溯源统计及 PDF |

[英文 README](../README.md)提供每一步可以单独复制执行的命令。

## 特征和模型

五种甲基化特征为 **MM、PDR、CHALM、MHL、Entropy**，另有 EM 和 BPM 两类片段组学特征。第一层每类特征分别训练 LR、GBM、线性 SVM、XGBoost 和 CatBoost；第二层用第一层的折外预测训练线性 SVM。癌症检测与组织溯源分别建模，均使用五折，完整示例会生成 360 个新的模型文件。

`--preset demo` 保留五算法、五折和两层结构，缩短内部调参与 boosting 迭代以便运行示例。`--preset study` 使用原脚本的训练设置。示例结果仅用于检查流程，不能作为论文的临床性能结果。

## marker 和自己的数据

示例 marker 位于 `examples/three_cancers/markers.tsv`，三列是 `task`、`channel`、`feature`。`task` 为 `detection` 或 `origin`；`channel` 为七类特征之一；甲基化 feature 为 `chr.start.end`，EM/BPM feature 为六碱基序列。每个任务、每个通道至少提供两个不重复的 marker，且必须存在于对应矩阵中。

模拟 BAM 用来展示单个样本的处理；三癌种模拟队列用于训练和独立验证。实际研究应对每个参与者的 BAM 分别执行 Step 1，再按事先确定的队列分组进入 Step 2。可以填写 `examples/cohorts_template.tsv` 并通过 `workflow/cohorts.py` 批量完成前两步。不同队列应使用相互独立的参与者。

验证和盲测的 `samples.tsv` 不含诊断标签；标签只存放在各自的 `truth.tsv` 中，Step 5 才会读取。示例文件对运行者可见，因此这里展示的是盲测的软件接口，实际研究的盲态管理还需独立执行。

运行结果保存在 `work/`，包括新训练的模型和各队列预测。仓库的 `.gitignore` 已排除该目录及 R 模型文件。随包附有依赖版本、输入说明、模拟 BAM 的预期输出和模拟队列的参考评估结果。
