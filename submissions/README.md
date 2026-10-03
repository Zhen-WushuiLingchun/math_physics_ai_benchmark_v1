# 匿名答卷与身份对应

最终总表采用下列 30 份 PDF；文件保留匿名编号，身份在评分完成后按参与者披露合并。跨题同字母只有在明确映射后才可视为同一配置。

| 模型配置 | 第 1 题 | 第 2 题 | 第 3 题 |
| --- | --- | --- | --- |
| GPT-6.1 Sol max | [s2_1b](supplement_2nd/01/s2_1b.pdf) | [s2_2b](supplement_2nd/02/s2_2b.pdf) | [s2_3b](supplement_2nd/03/s2_3b.pdf) |
| GPT-6.1 Sol xhigh | [s2_1a](supplement_2nd/01/s2_1a.pdf) | [s2_2a](supplement_2nd/02/s2_2a.pdf) | [s2_3a](supplement_2nd/03/s2_3a.pdf) |
| GPT-6 Astra high | [1b](original/01/1b.pdf) | [2e](original/02/2e.pdf) | [3e](original/03/3e.pdf) |
| GPT-6 Sol max | [s_1b](supplement/01/s_1b.pdf) | [2c](original/02/2c.pdf) | [3c](original/03/3c.pdf) |
| Opus 5.5 xhigh | [s_1a](supplement/01/s_1a.pdf) | [s_2a](supplement/02/s_2a.pdf) | [s_3a](supplement/03/s_3a.pdf) |
| GPT-6 Sol xhigh | [1c](original/01/1c.pdf) | [2d](original/02/2d.pdf) | [3d](original/03/3d.pdf) |
| GPT-5.6 Sol xhigh | [1a](original/01/1a.pdf) | [2a](original/02/2a.pdf) | [3a](original/03/3a.pdf) |
| GPT-6 Luna max | [s_1c](supplement/01/s_1c.pdf) | [2b](original/02/2b.pdf) | [3b](original/03/3b.pdf) |
| DeepSeek V4.1 Flash max | [s_1d](supplement/01/s_1d.pdf) | [s_2b](supplement/02/s_2b.pdf) | [s_3b](supplement/03/s_3b.pdf) |
| Grok 4.7 high | [1d](original/01/1d.pdf) | [2f](original/02/2f.pdf) | [3f](original/03/3f.pdf) |

original/ 保存原批次 19 份，supplement/ 保存第一轮补充 8 份，supplement_2nd/ 保存第二轮补充 6 份。Opus 原批次 1e、2g、3g 保留作追溯，最终采用 s_1a、s_2a、s_3a，原费用各计一次。

[PDF 清单与 SHA-256](../data/submissions.json)列出全部 33 份文件、采用状态及替代关系；[分项与费用 CSV](../data/scores.csv)包含 30 份采用结果。第二轮 a 为 GPT-6.1 Sol xhigh，b 为 GPT-6.1 Sol max，映射和费用来源见[成本记录](../costs/supplement_2nd_costs.json)。
