# 生物医学R包完整分类指南

## 🧬 **非编码RNA分析 | Non-coding RNA Analysis**
```
miRNAtap        - 多数据库miRNA靶基因预测 | Multiple database miRNA target prediction
multiMiR        - miRNA-mRNA互作数据库集成 | miRNA-mRNA interaction database integration
mirBaseConverter - miRNA ID转换工具 | miRNA ID conversion tool
lncRNAMap       - lncRNA功能注释数据库 | lncRNA functional annotation database
LncRNA2Target   - lncRNA靶基因预测 | lncRNA target gene prediction
LncRNADisease   - lncRNA疾病关联分析 | lncRNA disease association analysis
circRNAprofiler - circRNA识别和定量 | circRNA identification and quantification
circAtlas       - circRNA综合数据库 | Comprehensive circRNA database
tRNA            - tRNA序列和结构分析 | tRNA sequence and structure analysis
miRBaseVersions.db - miRNA版本数据库 | miRNA version database
```

## 🧪 **表观遗传学 | Epigenetics**
```
ChIPseeker      - ChIP-seq峰注释和可视化 | ChIP-seq peak annotation and visualization
DiffBind        - ChIP-seq差异结合分析 | ChIP-seq differential binding analysis
csaw            - ChIP-seq峰检测和计数 | ChIP-seq peak detection and counting
SICtools        - 表观遗传变异检测 | Epigenetic variation detection
IlluminaHumanMethylationEPICanno.ilmn10.hg4 - EPIC甲基化芯片注释 | EPIC methylation array annotation
IlluminaHumanMethylation450kanno.ilmn12.hg19 - 450K芯片注释 | 450K array annotation
EpiDISH         - 表观遗传细胞类型去卷积 | Epigenetic cell type deconvolution
MethylMix       - 甲基化驱动基因识别 | Methylation driver gene identification
minfi           - 甲基化微阵列分析 | Methylation microarray analysis
ChAMP           - 集成甲基化分析流程 | Integrated methylation analysis pipeline
ChAMPdata       - ChAMP示例数据 | ChAMP example data
wateRmelon      - 甲基化数据预处理 | Methylation data preprocessing
```

## 🩸 **免疫学分析 | Immunology Analysis**
```
immunarch       - 免疫组库数据分析和可视化 | Immune repertoire data analysis and visualization
VDJtools        - T细胞和B细胞受体分析 | T-cell and B-cell receptor analysis
Immunarch       - TCR/BCR谱系分析 | TCR/BCR repertoire analysis
ChangeO         - 免疫受体序列分析 | Immune receptor sequence analysis
shazam          - 体细胞超突变分析 | Somatic hypermutation analysis
Alakazam        - 免疫组库多样性分析 | Immune repertoire diversity analysis
immunedeconv    - 免疫细胞去卷积 | Immune cell deconvolution
ImmuLncRNA      - 免疫相关lncRNA分析 | Immune-related lncRNA analysis
Immunophenogram - 免疫表型分析 | Immunophenotype analysis
```

## 🔬 **单细胞多组学 | Single-cell Multi-omics**
```
scater          - 单细胞数据质量控制和可视化 | Single-cell data QC and visualization
scran           - 单细胞基因表达标准化 | Single-cell gene expression normalization
scDblFinder     - 双细胞检测 | Doublet detection
scds            - 双细胞识别算法 | Doublet identification algorithms
batchelor       - 单细胞批次校正 | Single-cell batch correction
zellkonverter   - 单细胞数据格式转换 | Single-cell data format conversion
scVI            - 单细胞变分推断 | Single-cell variational inference
scANVI          - 单细胞注释变分推断 | Single-cell annotated variational inference
totalVI         - 总RNA单细胞分析 | Total RNA single-cell analysis
Seurat          - 单细胞分析综合平台 | Comprehensive single-cell analysis platform
SeuratData      - Seurat示例数据集 | Seurat example datasets
SeuratDisk      - Seurat数据格式转换 | Seurat data format conversion
SeuratWrappers  - Seurat扩展功能 | Seurat extension functions
SingleCellExperiment - 单细胞实验数据容器 | Single-cell experiment data container
SingleR         - 单细胞自动注释 | Single-cell automatic annotation
slingshot       - 单细胞轨迹推断 | Single-cell trajectory inference
tradeSeq        - 轨迹差异表达 | Trajectory differential expression
monocle         - 单细胞轨迹分析 | Single-cell trajectory analysis
monocle3        - 新一代单细胞分析 | Next-generation single-cell analysis
AUCell          - 基因集活性评分 | Gene set activity scoring
CellChat        - 细胞间通信分析 | Cell-cell communication analysis
cicero          - 单细胞染色质可及性分析 | Single-cell chromatin accessibility analysis
```

## 🗺️ **空间转录组 | Spatial Transcriptomics**
```
SpatialExperiment - 空间转录组数据容器 | Spatial transcriptomics data container
STutility        - 空间转录组分析工具 | Spatial transcriptomics analysis tools
BayesSpace       - 贝叶斯空间聚类 | Bayesian spatial clustering
spotlight        - 细胞类型去卷积 | Cell type deconvolution
stLearn          - 空间转录组机器学习 | Spatial transcriptomics machine learning
SPOTlight        - 空间细胞类型映射 | Spatial cell type mapping
Giotto           - 空间转录组分析套件 | Spatial transcriptomics analysis suite
spacexr          - 空间表达数据分析 | Spatial expression data analysis
squidpy          - 空间分子分析Python接口 | Spatial molecular analysis Python interface
spatstat         - 空间点模式分析 | Spatial point pattern analysis
```

## 🤖 **深度学习与机器学习 | Deep Learning & Machine Learning**
```
keras           - R接口深度学习框架 | R interface to deep learning framework
tensorflow      - TensorFlow R接口 | TensorFlow R interface
torch           - PyTorch R接口 | PyTorch R interface
MXNet           - 高效深度学习框架 | Efficient deep learning framework
h2o             - 可扩展机器学习平台 | Scalable machine learning platform
caretEnsemble   - 集成机器学习模型 | Ensemble machine learning models
SuperLearner    - 超级学习器集成 | Super learner ensemble
mlr3            - 模块化机器学习框架 | Modular machine learning framework
mlr3learners    - mlr3学习器集合 | mlr3 learners collection
mlr3extralearners - 额外学习器扩展 | Extra learners extension
tidymodels      - 整洁建模框架 | Tidy modeling framework
parsnip         - 统一模型接口 | Unified model interface
workflows       - 建模工作流管理 | Modeling workflow management
tune            - 超参数调优 | Hyperparameter tuning
xgboost         - 梯度提升树 | Gradient boosting trees
lightgbm        - 轻量级梯度提升 | Lightweight gradient boosting
catboost        - 分类梯度提升 | Categorical gradient boosting
randomForest    - 随机森林 | Random forests
randomForestSRC - 随机森林生存分析 | Random forest survival analysis
ranger          - 快速随机森林 | Fast random forests
glmnet          - 正则化广义线性模型 | Regularized generalized linear models
gbm             - 梯度提升机 | Gradient boosting machine
e1071           - 支持向量机等算法 | Support vector machines and other algorithms
```

## 🧮 **生物统计学 | Biostatistics**
```
survival        - 生存分析基础 | Survival analysis foundation
survminer       - 生存分析可视化 | Survival analysis visualization
survivalROC     - 生存ROC分析 | Survival ROC analysis
survcomp        - 生存模型比较 | Survival model comparison
survMisc        - 生存分析工具 | Survival analysis tools
survRM2         - 限制平均生存时间 | Restricted mean survival time
survivalsvm     - 支持向量机生存分析 | Support vector machine survival analysis
rms             - 回归建模策略 | Regression modeling strategies
lme4            - 线性混合效应模型 | Linear mixed-effects models
nlme            - 非线性混合效应模型 | Nonlinear mixed-effects models
MCMCglmm        - MCMC广义线性混合模型 | MCMC generalized linear mixed models
brms            - 贝叶斯回归模型 | Bayesian regression models
rstanarm        - 斯坦回归模型 | Stan regression models
bayesplot       - 贝叶斯模型可视化 | Bayesian model visualization
multcomp        - 多重比较检验 | Multiple comparison tests
pwr             - 统计功效分析 | Statistical power analysis
pwr2            - 功效分析扩展 | Power analysis extension
faraway         - 实用回归分析 | Practical regression analysis
agricolae       - 农业实验统计 | Agricultural experimental statistics
DoE.base        - 实验设计基础 | Design of experiments basics
AlgDesign       - 实验设计算法 | Experimental design algorithms
```

## 📈 **时间序列分析 | Time Series Analysis**
```
zoo             - 规则和不规则时间序列 | Regular and irregular time series
xts             - 扩展时间序列处理 | Extended time series processing
tseries         - 时间序列分析和计算金融 | Time series analysis and computational finance
forecast        - 时间序列预测 | Time series forecasting
timetk          - 时间序列工具包 | Time series toolkit
```

## 🏥 **临床研究 | Clinical Research**
```
clinicalutils   - 临床数据分析工具 | Clinical data analysis tools
mediana         - 临床试验设计 | Clinical trial design
clinfun         - 临床相关函数 | Clinical functions
clinPK          - 临床药代动力学 | Clinical pharmacokinetics
clinReport      - 临床报告生成 | Clinical report generation
surveillance    - 疾病监测 | Disease surveillance
epiR            - 流行病学分析 | Epidemiological analysis
epitools        - 流行病学工具 | Epidemiology tools
incidence       - 发病率估计 | Incidence estimation
MatchIt         - 倾向得分匹配 | Propensity score matching
PredictABEL     - 疾病风险预测 | Disease risk prediction
riskRegression  - 风险回归模型 | Risk regression models
dynpred         - 动态预测模型 | Dynamic prediction models
pec             - 预测误差曲线 | Prediction error curves
ComparisonSurv  - 生存模型比较 | Survival model comparison
condsurv        - 条件生存分析 | Conditional survival analysis
```

## 🦠 **传染病与病毒学 | Infectious Diseases & Virology**
```
viral           - 病毒基因组分析 | Viral genome analysis
virusclust      - 病毒序列聚类 | Viral sequence clustering
phyloscan       - 系统发育扫描 | Phylogenetic scanning
OutbreakTools   - 疫情爆发数据分析 | Outbreak data analysis
EpiEstim        - 有效再生数估计 | Effective reproduction number estimation
EpiModel        - 传染病传播模型 | Infectious disease transmission modeling
Epi             - 流行病学数据分析 | Epidemiological data analysis
```

## 🧪 **蛋白质组学 | Proteomics**
```
MSnbase         - 质谱数据处理基础 | Mass spectrometry data processing foundation
MSstatsTMT      - TMT标记定量分析 | TMT labeling quantitative analysis
MSnID           - 质谱鉴定数据处理 | Mass spectrometry identification data processing
MSqRob          - 质谱定量稳健统计 | Mass spectrometry quantitative robust statistics
MZmine          - 质谱数据处理平台 | Mass spectrometry data processing platform
OpenMS          - 开源质谱系统 | Open source mass spectrometry system
proDA           - 蛋白质组学差异表达 | Proteomics differential expression
DEP             - 蛋白质组学分析流程 | Proteomics analysis pipeline
isobar          - 同位素标记定量 | Isotope labeling quantification
DIAassist       - DIA数据分析 | DIA data analysis
ProTraS         - 蛋白质翻译分析 | Protein translation analysis
protr           - 蛋白质序列特征计算 | Protein sequence feature computation
Rcintr          - 蛋白质相互作用中心 | Protein interaction centers
protViz         - 蛋白质组学可视化 | Proteomics visualization
Rdisop          - 蛋白质溶解度预测 | Protein solubility prediction
Bio3D           - 生物分子结构和序列 | Biomolecular structure and sequence
bio3d           - 蛋白质结构分析 | Protein structure analysis
```

## 🧬 **基因组学与序列分析 | Genomics & Sequence Analysis**
```
Biostrings      - 生物序列处理 | Biological sequence processing
IRanges         - 整数范围处理 | Integer range processing
GenomicRanges   - 基因组区间操作 | Genomic range operations
GenomicFeatures - 基因组特征注释 | Genomic feature annotation
XVector         - 高效序列存储 | Efficient sequence storage
S4Vectors       - 向量和列表数据结构 | Vector and list data structures
seqinr          - 生物序列分析 | Biological sequence analysis
msa             - 多序列比对 | Multiple sequence alignment
DNAcopy         - DNA拷贝数分析 | DNA copy number analysis
BSgenome        - 全基因组序列 | Whole genome sequences
BSgenome.Hsapiens.UCSC.hg19 - hg19参考基因组 | hg19 reference genome
BSgenome.Hsapiens.UCSC.hg38 - hg38参考基因组 | hg38 reference genome
ape             - 系统发育分析 | Phylogenetic analysis
phangorn        - 系统发育分析 | Phylogenetic analysis
geiger          - 进化分析 | Evolutionary analysis
phytools        - 系统发育工具 | Phylogenetic tools
phylobase       - 系统发育数据结构 | Phylogenetic data structures
```

## 🌐 **网络与系统生物学 | Network & Systems Biology**
```
igraph          - 网络分析和可视化 | Network analysis and visualization
igraphdata      - 图数据集 | Graph datasets
network         - 网络分析工具 | Network analysis tools
sna             - 社会网络分析 | Social network analysis
tnet            - 时间网络分析 | Temporal network analysis
multiplex       - 多层网络分析 | Multiplex network analysis
tidygraph       - 整洁图分析 | Tidy graph analysis
ggraph          - 图形语法网络可视化 | Grammar of graphics network visualization
networkD3       - D3网络可视化 | D3 network visualization
BioNet          - 生物网络分析 | Biological network analysis
NetPathMiner    - 网络路径挖掘 | Network path mining
Rgraph2D3       - R图转D3可视化 | R graph to D3 visualization
BoolNet         - 布尔网络分析 | Boolean network analysis
CellNOptR       - 细胞网络优化 | Cell network optimization
parcor          - 偏相关网络 | Partial correlation networks
CellDesigner    - 细胞通路建模 | Cell pathway modeling
synergyfinder   - 药物协同效应分析 | Drug synergy effect analysis
SynergyFinder   - 药物协同发现 | Drug synergy discovery
```

## 🧬 **功能富集与通路分析 | Functional Enrichment & Pathway Analysis**
```
clusterProfiler - 聚类基因功能富集 | Cluster gene functional enrichment
enrichplot      - 富集结果可视化 | Enrichment result visualization
DOSE            - 疾病本体富集 | Disease ontology enrichment
ReactomePA      - Reactome通路分析 | Reactome pathway analysis
pathview        - 通路数据可视化 | Pathway data visualization
KEGGgraph       - KEGG通路图处理 | KEGG pathway graph processing
KEGGREST        - KEGG数据库接口 | KEGG database interface
GOSemSim        - GO语义相似性 | GO semantic similarity
GO.db           - Gene Ontology数据库 | Gene Ontology database
fgsea           - 快速基因集富集分析 | Fast gene set enrichment analysis
GSEABase        - 基因集富集基础 | Gene set enrichment base
GSA             - 基因集分析 | Gene set analysis
SPIA            - 信号通路影响分析 | Signaling pathway impact analysis
```

## 📊 **数据可视化 | Data Visualization**
```
ggplot2         - 图形语法绘图系统 | Grammar of graphics plotting system
ggtree          - 系统树可视化 | Phylogenetic tree visualization
ggpubr          - 发表级图形 | Publication-ready graphics
ggsci           - 科学期刊配色 | Scientific journal color palettes
ggthemes        - 图形主题集合 | Graphics themes collection
ggtech          - 科技公司主题 | Tech company themes
gganimate       - 动画图形 | Animated graphics
rayshader       - 3D地形可视化 | 3D terrain visualization
rgl             - 3D交互图形 | 3D interactive graphics
plotrix         - 各种特殊图形 | Various special plots
ggforce         - ggplot2几何扩展 | ggplot2 geometry extensions
ggalt           - 替代坐标和统计 | Alternative coordinates and stats
ggExtra         - 边际图形 | Marginal graphics
ggformula       - 公式接口图形 | Formula interface graphics
ggstance        - 水平位置调整 | Horizontal position adjustments
ggrepel         - 智能标签防重叠 | Smart label repel
ggbio           - 基因组数据可视化 | Genomic data visualization
ComplexHeatmap  - 复杂热图绘制 | Complex heatmap drawing
pheatmap        - 热图可视化 | Heatmap visualization
plotly          - 交互式图形 | Interactive graphics
circlize        - 环形可视化 | Circular visualization
OmicCircos      - 组学环形图 | Omics circos plots
VennDiagram     - 韦恩图绘制 | Venn diagram drawing
ggVennDiagram   - ggplot2韦恩图 | ggplot2 Venn diagrams
ggalluvial      - 桑基图和冲积图 | Sankey and alluvial diagrams
gganatogram     - 解剖图可视化 | Anatogram visualization
ggbreak         - 坐标轴截断 | Axis breaking
ggcor           - 相关性图形 | Correlation graphics
ggcorrplot      - 相关矩阵可视化 | Correlation matrix visualization
ggnewscale      - 多标度图形 | Multiple scales graphics
ggpattern       - 图案填充 | Pattern fills
ggplotify       - 图形对象转换 | Graphics object conversion
ggpp            - 图形语法扩展 | Grammar of graphics extensions
ggprism         - Prism风格图形 | Prism-style graphics
cowplot         - 图形排列组合 | Plot arrangement and composition
patchwork       - 图形拼贴 | Plot patching
grid            - 网格图形系统 | Grid graphics system
gridBase        - 网格与传统图形集成 | Grid and base graphics integration
gridExtra       - 网格图形工具 | Grid graphics tools
RColorBrewer    - 调色板 | Color palettes
viridis         - 色盲友好配色 | Colorblind-friendly color scales
```

## 💾 **数据处理与操作 | Data Processing & Manipulation**
```
dplyr           - 数据操作语法 | Data manipulation grammar
tidyr           - 数据整理 | Data tidying
purrr           - 函数式编程 | Functional programming
readr           - 数据读取 | Data reading
readxl          - Excel文件读取 | Excel file reading
data.table      - 高效数据处理 | Efficient data processing
stringr         - 字符串处理 | String processing
janitor         - 数据清洗 | Data cleaning
reshape2        - 数据重塑 | Data reshaping
plyr            - 数据分割-应用-组合 | Split-apply-combine
disk.frame      - 磁盘大数据处理 | Disk-based big data processing
bigmemory       - 大内存矩阵 | Large memory matrices
ff              - 文件存储大数据 | File-based big data storage
bigstatsr       - 大规模统计计算 | Large-scale statistical computing
biglasso        - 大数据lasso回归 | Big data lasso regression
```

## 📋 **注释与元数据 | Annotation & Metadata**
```
AnnotationDbi   - 注释数据库接口 | Annotation database interface
AnnotationHub   - 注释资源中心 | Annotation resource hub
org.Hs.eg.db    - 人类基因注释 | Human gene annotation
org.Mm.eg.db    - 小鼠基因注释 | Mouse gene annotation
TxDb.Hsapiens.UCSC.hg19.knownGene - hg19转录本数据库 | hg19 transcript database
TxDb.Hsapiens.UCSC.hg38.knownGene - hg38转录本数据库 | hg38 transcript database
hgu133plus2.db  - 微阵列注释 | Microarray annotation
biomaRt         - BioMart数据库接口 | BioMart database interface
RIdeogram       - 染色体数据可视化 | Chromosome data visualization
```

## 🔧 **工作流与可重复研究 | Workflow & Reproducible Research**
```
rmarkdown       - 动态文档生成 | Dynamic document generation
knitr           - 报告生成引擎 | Report generation engine
bookdown        - 书籍创作 | Book authoring
blogdown        - 博客创建 | Blog creation
shiny           - 交互式Web应用 | Interactive web applications
flexdashboard   - 交互式仪表板 | Interactive dashboards
drake           - 数据流水线管理 | Data pipeline management
targets         - 流水线目标管理 | Pipeline target management
workflowr       - 可重复研究工作流 | Reproducible research workflows
renv            - 环境依赖性管理 | Environment dependency management
packrat         - 项目包管理 | Project package management
```

## 🧪 **药物发现与毒理学 | Drug Discovery & Toxicology**
```
pharmacoGx      - 药物基因组学 | Pharmacogenomics
PharmacoGx      - 药物敏感性分析 | Drug sensitivity analysis
xenograft       - 异种移植模型分析 | Xenograft model analysis
dr4pl           - 四参数剂量反应曲线 | Four-parameter dose-response curves
drc             - 剂量反应曲线分析 | Dose-response curve analysis
pRRophetic      - 药物反应预测 | Drug response prediction
pRRophetic2     - 改进的药物反应预测 | Improved drug response prediction
CMScaller       - 结直肠癌分子分型 | Colorectal cancer molecular subtyping
deconstructSigs - 突变特征分析 | Mutational signature analysis
MutationalPatterns - 突变模式分析 | Mutational patterns analysis
```

## 🔬 **高通量筛选 | High-throughput Screening**
```
cellHTS2        - 高通量筛选数据分析 | High-throughput screening data analysis
FlowSOM         - 流式数据自组织映射 | Flow data self-organizing maps
flowCore        - 流式数据基础结构 | Flow data infrastructure
flowWorkspace   - 流式数据分析环境 | Flow data analysis environment
CytoML          - 流式数据机器学习 | Flow data machine learning
openCyto        - 开放流式分析 | Open flow analysis
CytoExploreR    - 流式数据探索 | Flow data exploration
RchyOptimyx     - 流式细胞术优化 | Flow cytometry optimization
CATALYST        - 质谱流式数据分析 | Mass cytometry data analysis
flowAI          - 流式数据质量控制 | Flow data quality control
```

## 🌿 **微生物组与生态学 | Microbiome & Ecology**
```
phyloseq        - 微生物组数据分析 | Microbiome data analysis
microbiome      - 微生物组分析工具 | Microbiome analysis tools
vegan           - 群落生态学分析 | Community ecology analysis
metacoder       - 分类数据可视化 | Taxonomic data visualization
DEICODE         - 微生物组稀疏数据分析 | Microbiome sparse data analysis
ANCOMBC         - 微生物组差异丰度分析 | Microbiome differential abundance analysis
qurro           - 比例数据可视化 | Proportional data visualization
microbiomeMarker - 微生物组生物标志物 | Microbiome biomarkers
```

## 🧬 **癌症基因组学 | Cancer Genomics**
```
maftools        - 癌症突变数据分析 | Cancer mutation data analysis
TCGAbiolinks    - TCGA数据下载和分析 | TCGA data download and analysis
oncoPredict     - 药物敏感性预测 | Drug sensitivity prediction
CancerSubtypes  - 癌症亚型识别 | Cancer subtype identification
MuSiC           - 组织细胞类型分解 | Tissue cell type decomposition
EPIC            - 免疫细胞比例估计 | Immune cell proportion estimation
MCPcounter      - 免疫细胞计数 | Immune cell counting
ConsensusClusterPlus - 共识聚类 | Consensus clustering
iClusterPlus    - 多组学数据整合聚类 | Multi-omics integrative clustering
```

## 🛠️ **开发与工具包 | Development & Tools**
```
devtools        - 包开发工具 | Package development tools
usethis         - 自动化包设置 | Automated package setup
testthat        - 单元测试 | Unit testing
pkgdown         - 包文档网站 | Package documentation websites
roxygen2        - 文档生成 | Documentation generation
Rcpp            - R与C++集成 | R and C++ integration
reticulate      - Python接口 | Python interface
```

## 📦 **基础与核心包 | Foundation & Core Packages**
```
Biobase         - 生物实验数据基础 | Bioconductor experiment data base
BiocGenerics    - Bioconductor通用函数 | Bioconductor generic functions
BiocManager     - Bioconductor包管理 | Bioconductor package management
BiocParallel    - 并行计算 | Parallel computing
Matrix          - 稀疏矩阵处理 | Sparse matrix processing
utils           - 工具函数 | Utility functions
tools           - 包开发工具 | Package development tools
stats           - 统计函数 | Statistical functions
graphics        - 基础图形 | Base graphics
grDevices       - 图形设备 | Graphics devices
```
