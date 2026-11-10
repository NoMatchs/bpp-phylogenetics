# Bio++ Phylogenetics Analysis

这是我的 Bio++ 学习与实践项目。记录了从系统安装、C++ 开发到最大似然建树的完整流程。

## 📊 数据来源
* 数据从 TreeBASE 数据库获取（Matrix M16588）。
* 包含 10 种哺乳动物，共 590 个碱基的比对数据（内含子序列）。
* 在分析前，通过 Bio++ 提取了 163 个完全位点（不含 Gap）。

## 🛠️ 环境与依赖
* Ubuntu 24.04
* 从源码编译安装了 Bio++ 3.0 (Core + Seq + Phyl)
* Eigen3, Biopython

## 🚀 功能与代码结构
* `real_alignment.fasta`: 原始多序列比对数据。
* `ml_real.bpp`: bppml 命令行参数文件（GTR + Gamma 模型）。
* `ml_cpp_test.cpp`: 使用 Bio++ 3.0 全新 DataFlow 架构计算的 C++ 测试代码。
* `ml_real_output.tree`: 最大似然法优化后的 Newick 格式进化树。

## 📈 分析结果
在 GTR+Gamma 模型下，估计出的 Gamma 形状参数 `alpha = 1.567`，证明了序列中存在显著的位点间速率异质性。

最终的无根树（可通过 FigTree 查看）展示了哺乳动物的系统发育关系，其中食肉目（狗、猫）与鲸偶蹄目（海豚、牛、猪）的分支结构非常清晰。
