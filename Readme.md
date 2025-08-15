## 简介

从`fasta`格式的核酸序列开始，执行 VDJ gene annotation、CDR detection、SHM count，以获取 B 细胞受体（BCR）的各种信息。使用`igblast, changeo, ANARCI`三个软件互相矫正，获取最终结果

## 示例数据

### 准备

首先为系统安装 conda，推荐使用 [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) 发行版

使用[chsrc](https://github.com/RubyMetric/chsrc)，将 conda 源切换到最合适的镜像站，尽量避免环境安装失败

重新连接终端后，安装 Snakemake 基础 conda 环境

```
mamba create -n vdj && mamba activate vdj
mamba install snakemake=8.24.0 just=1.42.4 -y
```

### 开始运行

```
snakemake --profile profiles/local --configfile configs/eg.yaml -j1
```

该命令将自动完成所需匿名 conda 环境的安装，并运行示例数据

## 特性

使用`igblast`运行结果，判断序列为轻链或重链

使用`configfile`中`cell_pattern`，正则匹配提取轻重链序列对应的细胞。默认`'^(.+?)-[HLK]$'`，例如`.fasta`文件中，`>cell1-H`和`>cell1-L`被认定为 cell1，`>cell9-H`和`>cell9-K`被认定为 cell9

与 `cell10x` pipeline 相比：

- 使用 igblast 的`sequence`经过 biostrings 翻译，获取氨基酸序列，然后跑 ANARCI

- `seq_aa`从 igblast 获取，`seq_align_aa`仍然从 ANARCI 获取

测试集：

- `eg1` 是常规测试集

- `eg2` 是 igblast 的`seq_aa`无法被 ANARCI 处理的代表情况。但经 biostrings 翻译后，可修复该问题
