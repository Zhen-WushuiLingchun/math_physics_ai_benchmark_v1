# 汇总与复现脚本

从仓库根目录运行。Python 3.10+；绘图依赖列在 `requirements.txt`。

```bash
python scripts/check_dataset.py
python -m pip install -r scripts/requirements.txt
python scripts/plots.py --conclusion-only
python scripts/write_report.py
```

`check_dataset.py` 使用标准库检查 27 份答卷哈希、24 项评分、CSV、三题总分与成本、前沿和工作簿源单元格。它不会重新裁定数学物理结论。

`plots.py --conclusion-only` 绘制最终对数／线性结论图；不带参数还绘制逐题图、预算组合图，并读取现有工作簿核对图表数据。需安装支持中文的字体，可设置 `BENCHMARK_FONT` 为字体文件路径。所有图片写入 `results/figures/`。

`write_report.py` 用完整数据生成 `results/analysis.md`，模型评价文字来自既有评阅。`build_workbook.mjs` 使用 `@oai/artifact-tool`，需要提供该包的 Codex artifact runtime；它不是仅用上述 Python 依赖即可运行的脚本。具备该运行环境时：

```bash
node scripts/build_workbook.mjs analyze
node scripts/build_workbook.mjs author
```

工作簿作者模式在现有 Excel 文件上更新三张 Task 表和综合汇总。统计公式关联 Task 表；原生散点图的横坐标为生成时的成本快照，费用变动后应重新生成图表。暂存核验输出写入被 Git 忽略的 `.work/`。

`evaluations/` 和 `runs/` 内的科学计算属于评阅者或作答模型的独立材料，不会由这些汇总脚本自动重跑。
