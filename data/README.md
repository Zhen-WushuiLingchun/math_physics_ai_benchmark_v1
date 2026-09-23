# 数据口径与文件

- `analysis_data.json`：当前完整数据，`rows` 为采用的 24 份回答，`full` 为三题等权汇总，包含逐题前沿、512 个组合中的有效前沿及代表方案。JSON 中文件路径相对仓库根目录。
- `scores.csv`：同一批 24 项的扁平表，含匿名编号、模型、五维分数、总分、评阅范围、参考成本、PDF 和替代关系。
- `submissions.json`：全部 27 份 PDF 的 SHA-256 与采用状态。
- `original_model_grading.xlsx`：收到的原始模型—费用表，保留作输入记录；最终身份纠正、补充评分和汇总以 `analysis_data.json` 与 `results/` 中的工作簿为准。
- `input_rows.json`：此前工作簿原有费用记录，用于重建原批次费用和 Opus 费用复用。
- `chart_bindings.json`：工作簿中四个原生散点图的坐标与单元格对应。
- `qa_checks.json`：工作簿导出检查记录。

分数按既有评阅结果保留；成本沿用提供者记录，单位 USD。各题相同字母不代表同一个模型。评阅范围不是统计置信区间。评分和费用信息应与根目录 README 的方法限制一起使用。
