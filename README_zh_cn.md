# GeneMiner2 UCE fork: 面向UCE流程的GeneMiner2扩展分支
**[View README in English](README.md)**

# 介绍
本仓库是面向target-enrichment和UCE数据流程的GeneMiner2命令行扩展分支。该分支保留原版GeneMiner2基于参考序列的分子标记恢复流程，同时增加了更适合UCE侧翼序列延伸、phyluce衔接和比对过滤的命令行选项。

本分支的主要新增内容包括：
- `--assembly-mode uce`：放宽参考边界裁剪，并优先保留更长且有reads支持的UCE侧翼序列。
- 在`uce_contigs/`下导出phyluce兼容的UCE contig文件。
- 通过`--alignment-filter alifilter`支持可选的AliFilter比对列过滤。
- 通过`--msa-threads`和`--filter-processes`改进combine阶段的并行控制。

原版GeneMiner2是一款为系统发育基因组学设计的全功能工具包，软件主要功能包括：
- 从二代测序数据中挖掘单拷贝核基因、质体基因等分子标记
- 将多个分子标记切齐、排序、建立串联和溯祖系统发育树
- 拼接注释动植物质体基因组

用户能够在GeneMiner2中完成从NGS数据获取到系统发育树建立的所有工作。本分支的UCE相关改动目前只针对命令行流程，图形界面未修改。

![](images/fig1.png)

![](images/fig2.jpg)

## 引用

欢迎使用并引用GeneMiner2:

Yu XY, Tang ZZ, Zhang Z, Song YX, He H, Shi Y, Hou JQ, Yu Y. 2026. **GeneMiner2**: Accurate and automated recovery of genes from genome-skimming data. *Molecular Ecology Resources* 26: e70111.https://doi.org/10.1111/1755-0998.70111



相关前期分子标记挖掘工具，欢迎引用:

Zhang Z, Xie PL, Guo YL, Zhou WB, Liu EY, Yu Y. 2022. **Easy353**: A tool to get Angiosperms353 genes for phylogenomic research. *Molecular Biology and Evolution* 39(12): msac261.https://doi.org/10.1093/molbev/msac261


Xie PL, Guo YL, Teng Y, Zhou WB, Yu Y. 2024. **GeneMiner**: A tool for extracting phylogenetic markers from next-generation sequencing data. *Molecular Ecology Resources* 24(3): e13924.https://doi.org/10.1111/1755-0998.13924

如果使用`--alignment-filter alifilter`，也请引用：

Bianchini G, Zhu R, Cicconardi F, Moody ERR. 2026. **AliFilter: a machine learning approach to alignment filtering.** *Molecular Biology and Evolution* 43(4): msag097. https://doi.org/10.1093/molbev/msag097


# 安装和需求

原版GeneMiner2安装包可从下面的地址获取:

**[SourceForge](https://sourceforge.net/projects/geneminer/files/)**

SourceForge上的安装包属于原版GeneMiner2发布版本，不一定包含本分支新增的UCE、AliFilter和combine阶段性能参数。如需使用本分支的命令行功能，请克隆本仓库并直接运行脚本或从源码构建CLI。

## Windows用户

请从上面的Sourceforge地址下载最新的Windows软件包（**GeneMiner_win64_XXXXXXXX.zip**）并解压到一个独立的文件夹中。双击GeneMiner.exe启动图形界面。

**注意:** 请勿在移动磁盘上运行GeneMiner2。请勿将GeneMiner2安装到带有中文的文件夹中。

**关于并行运行:** 不要在同一个文件夹下打开多个窗口并行运行，可以复制GeneMiner所在的文件夹，并在副本文件中运行第二个窗口。文件夹路径不能有中文名。

## macOS用户

请从上面的Sourceforge地址下载最新的macOS安装映像（**GeneMiner_macos_XXXXXXXX.dmg**），并将GeneMiner.app拖放到想要的位置。

如果遇到 "GeneMiner.app已损坏，无法打开"之类的错误，请打开终端，运行类似下面的命令:
- xattr -cr **把GeneMiner.app拖放到此处**
- 例如: xattr -cr /Applications/GeneMiner.app

由于GeneMiner2通过Wine技术运行在macOS上，运行效率远低于Windows，因此不建议用GeneMiner2在macOS上进行大规模的分析。

## Linux和其他\*nix用户

在桌面Linux环境下，请考虑使用兼容性工具运行Windows版本：

**[在Linux上运行Windows版本](manual/ZH_CN/linux_desktop.md)**

原版GeneMiner2的服务器Linux用户可以从上面的Sourceforge地址下载Linux命令行版本（**GeneMiner_cli_linux_XXXXXXXX.tar.gz**）并解压。这一版本支持Debian 11或以上、Ubuntu 20.04或以上、AlmaLinux 9或以上。此外，请确保你安装了libbz2、libgomp和zlib依赖库。在Ubuntu上，可以用这行命令安装这些依赖：

```bash
sudo apt-get install libbz2 libgomp1 zlib1g
```

**[命令行版本的使用说明](manual/ZH_CN/command_line.md#用法)**

如需使用本分支的UCE、AliFilter和combine阶段性能参数，请从本仓库手动编译命令行版本，或直接使用`scripts`文件夹中的Python脚本。这些脚本提供了GeneMiner2的核心命令行功能，可以在任何操作系统上部署。

**[从头编译命令行版本](manual/ZH_CN/command_line.md)**

---

# 示例演示

[Tutorial 1 - 快速掌握使用方法(单个样品)](/DEMO/DEMO1/DEMO1.md)


[Tutorial 2 - 获取质体基因组和质体基因](/DEMO/DEMO2/DEMO2.md)

  
[Tutorial 3 - 单拷贝基因建树流程(批量，多个样品)](DEMO/DEMO3/DEMO3.md)





# 详细说明


原版完整使用手册请移步 [manual](manual/manual_geneminer.pdf)。本分支新增的命令行选项以[命令行版本的使用说明](manual/ZH_CN/command_line.md)为准。

有关软件图形界面和功能的详细说明[请见此处](manual/ZH_CN/readme_detailed.md)

有关输出目录详解[请见此处](manual/ZH_CN/output.md)


# 常见问题
[常见问题解答](manual/ZH_CN/FAQ.md)


# 联系方式
有关GeneMiner任何建议、问题，请联系邮箱
Xinyi_Yu2021@163.com.
