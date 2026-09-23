# 模型作答工程

本目录保留此前公开的模型作答工程，包括报告源码、计算代码和小型输出，目录名沿用原实验名称，内部相对布局保持不变。模型答卷包含的科学结论不因收录而视为已验证，正确性和适用范围请对照[评阅记录](../evaluations/README.md)。

| 配置 | 已有可编辑工程 |
| --- | --- |
| GPT-6 Astra high | [题 1](gpt_test_6_astra/solution/README.md) |
| GPT-6 Sol max | [题 1](Sol_max_test/) |
| GPT-6 Sol xhigh | [题 1–3](gpt_6_sol/) |
| GPT-6 Luna max | [题 1](luna_test/) |
| GPT-5.6 Sol xhigh | [题 1–3](gpt_test2_5.6sol/) |
| Opus 5.5 xhigh | [题 1](opus5.5_test/q1/) · [题 2](opus5.5_test_2/) · [题 3](opus5.5_test_3/) |
| DeepSeek V4.1 Flash max | [题 1](ds_flash_test/test1/) · [题 2](ds_flash_test2/) · [题 3](ds_flash_test3/) |
| Grok 4.7 high | [题 1–3](grok_test_4.7/) |

上述源码覆盖范围与最终答卷 PDF 的覆盖范围不同。所有 8 个配置的三题采用 PDF 均在 [submissions](../submissions/README.md)，未随本次资料提供的源码不作补写。

例如，从仓库根目录运行 Astra 第 1 题报告附带的精确检查：

```bash
python runs/gpt_test_6_astra/solution/checks/verify_q1.py
```

需要 Python 和 SymPy，详见该工程的 [README](gpt_test_6_astra/solution/README.md)。其他工程的 Python、Wolfram Language 和 XeLaTeX 依赖各不相同。目录整理没有重新执行全部科学计算，也不改变匿名评分时实际可获得的材料。
