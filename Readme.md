与 cell10x 相比：

- 使用 igblast 的`sequence`经过 biostrings 翻译，获取氨基酸序列，然后跑 ANARCI

- `seq_aa`从 igblast 获取，`seq_align_aa`仍然从 ANARCI 获取

测试集：

- `eg1` 是常规测试集

- `eg2` 是 igblast 的`seq_aa`无法被 ANARCI 处理的代表情况。但经 biostrings 翻译后，可修复该问题
