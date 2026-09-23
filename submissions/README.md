# 匿名答卷与身份对应

最终总表采用下列 24 份 PDF；文件保留匿名编号，身份在评分完成后按参与者提供的信息合并。每一列中的相同字母不代表跨题身份一致。

| 模型配置 | 第 1 题 | 第 2 题 | 第 3 题 |
| --- | --- | --- | --- |
| GPT-6 Astra high | [1b](original/01/1b.pdf) | [2e](original/02/2e.pdf) | [3e](original/03/3e.pdf) |
| GPT-6 Sol max | [s_1b](supplement/01/s_1b.pdf) | [2c](original/02/2c.pdf) | [3c](original/03/3c.pdf) |
| Opus 5.5 xhigh | [s_1a](supplement/01/s_1a.pdf) | [s_2a](supplement/02/s_2a.pdf) | [s_3a](supplement/03/s_3a.pdf) |
| GPT-6 Sol xhigh | [1c](original/01/1c.pdf) | [2d](original/02/2d.pdf) | [3d](original/03/3d.pdf) |
| GPT-5.6 Sol xhigh | [1a](original/01/1a.pdf) | [2a](original/02/2a.pdf) | [3a](original/03/3a.pdf) |
| GPT-6 Luna max | [s_1c](supplement/01/s_1c.pdf) | [2b](original/02/2b.pdf) | [3b](original/03/3b.pdf) |
| DeepSeek V4.1 Flash max | [s_1d](supplement/01/s_1d.pdf) | [s_2b](supplement/02/s_2b.pdf) | [s_3b](supplement/03/s_3b.pdf) |
| Grok 4.7 high | [1d](original/01/1d.pdf) | [2f](original/02/2f.pdf) | [3f](original/03/3f.pdf) |

`original/` 保存原批次 19 份，`supplement/` 保存补充批次 8 份。Opus 的原批次 `1e`、`2g`、`3g` 保留作追溯，最终采用对应的 `s_1a`、`s_2a`、`s_3a`，费用沿用原记录且各计一次。

[PDF 清单与 SHA-256](../data/submissions.json)列出全部 27 份文件、模型、题号、采用状态及替代关系。[最终分项和费用 CSV](../data/scores.csv)只包含采用的 24 份结果。
