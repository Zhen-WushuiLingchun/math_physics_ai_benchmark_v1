# 评分与核验记录

评分由食司使用 GPT 6 Astra Max 完成。题 1 由主评处理，题 2、3 使用 Astra max 子智能体评阅并由主评复核汇总与关键结论。评分依据是[同一份五维细则](../questions/math_physics_ai_benchmark_rubric.md)。

| 批次 | 第 1 题 | 第 2 题 | 第 3 题 | 机读汇总 |
| --- | --- | --- | --- | --- |
| 原批次 | [评阅](original/q1_assessment.md) | [评阅](original/q2_assessment.md) | [评阅](original/q3_assessment.md) | [19 份记录](original/评分明细.json) |
| 补充批次 | [评阅](supplement/q1_assessment.md) | [评阅](supplement/q2_assessment.md) | [评阅](supplement/q3_assessment.md) | [8 份记录](supplement/补充评分明细.json) |

批次报告保留评分当时的叙述、排名和匿名阶段记录，不单独作为当前完整总榜。当前采用版本、模型身份与成本见[最终数据](../data/analysis_data.json)和[完整分析](../results/analysis.md)。三份 Opus 补充答卷用于最终统计，原批次仍保留以便追溯。

评分时向模型隐去作答型号，不排除同一族模型在出题和评分时的隐性偏好。补充批次记录注明辅助评阅接触 `s_3a` 潜在 PDF 元数据身份线索的情况，主评未接触该线索并独立复核其主论证。评阅范围是合理复评范围，不是统计置信区间。

评阅者计算保存在 `original/q2/`、`original/q3/`、`supplement/q3_work/` 及 [q1_verification_record.json](supplement/q1_verification_record.json)；这些证据不计为作答模型自行完成的验证。保留的数值结果不是所有符号工具会话的完整执行日志。

`supplement/assemble_supplement.py` 是匿名批次汇总的历史脚本，会初始化该批次身份与成本字段；最终发布数据的入口是 `data/analysis_data.json`，不要用历史脚本覆盖已经整合身份的数据。
