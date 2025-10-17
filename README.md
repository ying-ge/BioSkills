# 生物医学R包完整分类指南

## 🧬 **非编码RNA分析**
```
miRNAtap        - 多数据库miRNA靶基因预测
multiMiR        - miRNA-mRNA互作数据库集成
mirBaseConverter - miRNA ID转换工具
lncRNAMap       - lncRNA功能注释数据库
LncRNA2Target   - lncRNA靶基因预测
LncRNADisease   - lncRNA疾病关联分析
circRNAprofiler - circRNA识别和定量
circAtlas       - circRNA综合数据库
tRNA            - tRNA序列和结构分析
miRBaseVersions.db - miRNA版本数据库
```

## 🧪 **表观遗传学**
```
ChIPseeker      - ChIP-seq峰注释和可视化
DiffBind        - ChIP-seq差异结合分析
csaw            - ChIP-seq峰检测和计数
SICtools        - 表观遗传变异检测
IlluminaHumanMethylationEPICanno.ilmn10.hg4 - EPIC甲基化芯片注释
IlluminaHumanMethylation450kanno.ilmn12.hg19 - 450K芯片注释
EpiDISH         - 表观遗传细胞类型去卷积
MethylMix       - 甲基化驱动基因识别
minfi           - 甲基化微阵列分析
ChAMP           - 集成甲基化分析流程
ChAMPdata       - ChAMP示例数据
wateRmelon      - 甲基化数据预处理
```

## 🩸 **免疫学分析**
```
immunarch       - 免疫组库数据分析和可视化
VDJtools        - T细胞和B细胞受体分析
Immunarch       - TCR/BCR谱系分析
ChangeO         - 免疫受体序列分析
shazam          - 体细胞超突变分析
Alakazam        - 免疫组库多样性分析
immunedeconv    - 免疫细胞去卷积
ImmuLncRNA      - 免疫相关lncRNA分析
Immunophenogram - 免疫表型分析
```

## 🔬 **单细胞多组学**
```
scater          - 单细胞数据质量控制和可视化
scran           - 单细胞基因表达标准化
scDblFinder     - 双细胞检测
scds            - 双细胞识别算法
batchelor       - 单细胞批次校正
zellkonverter   - 单细胞数据格式转换
scVI            - 单细胞变分推断
scANVI          - 单细胞注释变分推断
totalVI         - 总RNA单细胞分析
Seurat          - 单细胞分析综合平台
SeuratData      - Seurat示例数据集
SeuratDisk      - Seurat数据格式转换
SeuratWrappers  - Seurat扩展功能
SingleCellExperiment - 单细胞实验数据容器
SingleR         - 单细胞自动注释
slingshot       - 单细胞轨迹推断
tradeSeq        - 轨迹差异表达
monocle         - 单细胞轨迹分析
monocle3        - 新一代单细胞分析
AUCell          - 基因集活性评分
CellChat        - 细胞间通信分析
cicero          - 单细胞染色质可及性分析
```

## 🗺️ **空间转录组**
```
SpatialExperiment - 空间转录组数据容器
STutility        - 空间转录组分析工具
BayesSpace       - 贝叶斯空间聚类
spotlight        - 细胞类型去卷积
stLearn          - 空间转录组机器学习
SPOTlight        - 空间细胞类型映射
Giotto           - 空间转录组分析套件
spacexr          - 空间表达数据分析
squidpy          - 空间分子分析Python接口
spatstat         - 空间点模式分析
```

## 🤖 **深度学习与机器学习**
```
keras           - R接口深度学习框架
tensorflow      - TensorFlow R接口
torch           - PyTorch R接口
MXNet           - 高效深度学习框架
h2o             - 可扩展机器学习平台
caretEnsemble   - 集成机器学习模型
SuperLearner    - 超级学习器集成
mlr3            - 模块化机器学习框架
mlr3learners    - mlr3学习器集合
mlr3extralearners - 额外学习器扩展
tidymodels      - 整洁建模框架
parsnip         - 统一模型接口
workflows       - 建模工作流管理
tune            - 超参数调优
xgboost         - 梯度提升树
lightgbm        - 轻量级梯度提升
catboost        - 分类梯度提升
randomForest    - 随机森林
randomForestSRC - 随机森林生存分析
ranger          - 快速随机森林
glmnet          - 正则化广义线性模型
gbm             - 梯度提升机
e1071           - 支持向量机等算法
```

## 🧮 **生物统计学**
```
survival        - 生存分析基础
survminer       - 生存分析可视化
survminer       - 生存曲线和风险表
survivalROC     - 生存ROC分析
survcomp        - 生存模型比较
survMisc        - 生存分析工具
survRM2         - 限制平均生存时间
survivalsvm     - 支持向量机生存分析
rms             - 回归建模策略
lme4            - 线性混合效应模型
nlme            - 非线性混合效应模型
MCMCglmm        - MCMC广义线性混合模型
brms            - 贝叶斯回归模型
rstanarm        - 斯坦回归模型
bayesplot       - 贝叶斯模型可视化
multcomp        - 多重比较检验
pwr             - 统计功效分析
pwr2            - 功效分析扩展
faraway         - 实用回归分析
agricolae       - 农业实验统计
DoE.base        - 实验设计基础
AlgDesign       - 实验设计算法
```

## 📈 **时间序列分析**
```
zoo             - 规则和不规则时间序列
xts             - 扩展时间序列处理
tseries         - 时间序列分析和计算金融
forecast        - 时间序列预测
timetk          - 时间序列工具包
```

## 🏥 **临床研究**
```
clinicalutils   - 临床数据分析工具
mediana         - 临床试验设计
clinfun         - 临床相关函数
clinPK          - 临床药代动力学
clinReport      - 临床报告生成
surveillance    - 疾病监测
epiR            - 流行病学分析
epitools        - 流行病学工具
incidence       - 发病率估计
MatchIt         - 倾向得分匹配
rms             - 临床预测模型
PredictABEL     - 疾病风险预测
riskRegression  - 风险回归模型
dynpred         - 动态预测模型
pec             - 预测误差曲线
ComparisonSurv  - 生存模型比较
condsurv        - 条件生存分析
```

## 🦠 **传染病与病毒学**
```
viral           - 病毒基因组分析
virusclust      - 病毒序列聚类
phyloscan       - 系统发育扫描
OutbreakTools   - 疫情爆发数据分析
EpiEstim        - 有效再生数估计
EpiModel        - 传染病传播模型
Epi             - 流行病学数据分析
```

## 🧪 **蛋白质组学**
```
MSnbase         - 质谱数据处理基础
MSstatsTMT      - TMT标记定量分析
MSnID           - 质谱鉴定数据处理
MSqRob          - 质谱定量稳健统计
MZmine          - 质谱数据处理平台
OpenMS          - 开源质谱系统
proDA           - 蛋白质组学差异表达
DEP             - 蛋白质组学分析流程
isobar          - 同位素标记定量
DIAassist       - DIA数据分析
ProTraS         - 蛋白质翻译分析
protr           - 蛋白质序列特征计算
Rcintr          - 蛋白质相互作用中心
protViz         - 蛋白质组学可视化
Rdisop          - 蛋白质溶解度预测
Bio3D           - 生物分子结构和序列
bio3d           - 蛋白质结构分析
```

## 🧬 **基因组学与序列分析**
```
Biostrings      - 生物序列处理
IRanges         - 整数范围处理
GenomicRanges   - 基因组区间操作
GenomicFeatures - 基因组特征注释
XVector         - 高效序列存储
S4Vectors       - 向量和列表数据结构
seqinr          - 生物序列分析
msa             - 多序列比对
DNAcopy         - DNA拷贝数分析
BSgenome        - 全基因组序列
BSgenome.Hsapiens.UCSC.hg19 - hg19参考基因组
BSgenome.Hsapiens.UCSC.hg38 - hg38参考基因组
ape             - 系统发育分析
phangorn        - 系统发育分析
geiger          - 进化分析
phytools        - 系统发育工具
phylobase       - 系统发育数据结构
```

## 🌐 **网络与系统生物学**
```
igraph          - 网络分析和可视化
igraphdata      - 图数据集
network         - 网络分析工具
sna             - 社会网络分析
tnet            - 时间网络分析
multiplex       - 多层网络分析
tidygraph       - 整洁图分析
ggraph          - 图形语法网络可视化
networkD3       - D3网络可视化
BioNet          - 生物网络分析
NetPathMiner    - 网络路径挖掘
Rgraph2D3       - R图转D3可视化
BoolNet         - 布尔网络分析
CellNOptR       - 细胞网络优化
parcor          - 偏相关网络
CellDesigner    - 细胞通路建模
synergyfinder   - 药物协同效应分析
SynergyFinder   - 药物协同发现
```

## 🧬 **功能富集与通路分析**
```
clusterProfiler - 聚类基因功能富集
enrichplot      - 富集结果可视化
DOSE            - 疾病本体富集
ReactomePA      - Reactome通路分析
pathview        - 通路数据可视化
KEGGgraph       - KEGG通路图处理
KEGGREST        - KEGG数据库接口
GOSemSim        - GO语义相似性
GO.db           - Gene Ontology数据库
fgsea           - 快速基因集富集分析
GSEABase        - 基因集富集基础
GSA             - 基因集分析
SPIA            - 信号通路影响分析
```

## 📊 **数据可视化**
```
ggplot2         - 图形语法绘图系统
ggtree          - 系统树可视化
ggpubr          - 发表级图形
ggsci           - 科学期刊配色
ggthemes        - 图形主题集合
ggtech          - 科技公司主题
gganimate       - 动画图形
rayshader       - 3D地形可视化
rgl             - 3D交互图形
plotrix         - 各种特殊图形
ggforce         - ggplot2几何扩展
ggalt           - 替代坐标和统计
ggExtra         - 边际图形
ggformula       - 公式接口图形
ggstance        - 水平位置调整
ggrepel         - 智能标签防重叠
ggbio           - 基因组数据可视化
ComplexHeatmap  - 复杂热图绘制
pheatmap        - 热图可视化
plotly          - 交互式图形
circlize        - 环形可视化
OmicCircos      - 组学环形图
VennDiagram     - 韦恩图绘制
ggVennDiagram   - ggplot2韦恩图
ggalluvial      - 桑基图和冲积图
gganatogram     - 解剖图可视化
ggbreak         - 坐标轴截断
ggcor           - 相关性图形
ggcorrplot      - 相关矩阵可视化
ggnewscale      - 多标度图形
ggpattern       - 图案填充
ggplotify       - 图形对象转换
ggpp            - 图形语法扩展
ggprism         - Prism风格图形
cowplot         - 图形排列组合
patchwork       - 图形拼贴
grid            - 网格图形系统
gridBase        - 网格与传统图形集成
gridExtra       - 网格图形工具
RColorBrewer    - 调色板
viridis         - 色盲友好配色
```

## 💾 **数据处理与操作**
```
dplyr           - 数据操作语法
tidyr           - 数据整理
purrr           - 函数式编程
readr           - 数据读取
readxl          - Excel文件读取
data.table      - 高效数据处理
stringr         - 字符串处理
janitor         - 数据清洗
reshape2        - 数据重塑
plyr            - 数据分割-应用-组合
disk.frame      - 磁盘大数据处理
bigmemory       - 大内存矩阵
ff              - 文件存储大数据
bigstatsr       - 大规模统计计算
biglasso        - 大数据lasso回归
```

## 📋 **注释与元数据**
```
AnnotationDbi   - 注释数据库接口
AnnotationHub   - 注释资源中心
org.Hs.eg.db    - 人类基因注释
org.Mm.eg.db    - 小鼠基因注释
TxDb.Hsapiens.UCSC.hg19.knownGene - hg19转录本数据库
TxDb.Hsapiens.UCSC.hg38.knownGene - hg38转录本数据库
hgu133plus2.db  - 微阵列注释
biomaRt         - BioMart数据库接口
RIdeogram       - 染色体数据可视化
```

## 🔧 **工作流与可重复研究**
```
rmarkdown       - 动态文档生成
knitr           - 报告生成引擎
bookdown        - 书籍创作
blogdown        - 博客创建
shiny           - 交互式Web应用
flexdashboard   - 交互式仪表板
drake           - 数据流水线管理
targets         - 流水线目标管理
workflowr       - 可重复研究工作流
renv            - 环境依赖性管理
packrat         - 项目包管理
```

## 🧪 **药物发现与毒理学**
```
pharmacoGx      - 药物基因组学
PharmacoGx      - 药物敏感性分析
xenograft       - 异种移植模型分析
dr4pl           - 四参数剂量反应曲线
drc             - 剂量反应曲线分析
pRRophetic      - 药物反应预测
pRRophetic2     - 改进的药物反应预测
CMScaller       - 结直肠癌分子分型
deconstructSigs - 突变特征分析
MutationalPatterns - 突变模式分析
```

## 🔬 **高通量筛选**
```
cellHTS2        - 高通量筛选数据分析
FlowSOM         - 流式数据自组织映射
flowCore        - 流式数据基础结构
flowWorkspace   - 流式数据分析环境
CytoML          - 流式数据机器学习
openCyto        - 开放流式分析
CytoExploreR    - 流式数据探索
RchyOptimyx     - 流式细胞术优化
CATALYST        - 质谱流式数据分析
flowAI          - 流式数据质量控制
```

## 🌿 **微生物组与生态学**
```
phyloseq        - 微生物组数据分析
microbiome      - 微生物组分析工具
vegan           - 群落生态学分析
metacoder       - 分类数据可视化
DEICODE         - 微生物组稀疏数据分析
ANCOMBC         - 微生物组差异丰度分析
qurro           - 比例数据可视化
microbiomeMarker - 微生物组生物标志物
```

## 🧬 **癌症基因组学**
```
maftools        - 癌症突变数据分析
TCGAbiolinks    - TCGA数据下载和分析
oncoPredict     - 药物敏感性预测
CancerSubtypes  - 癌症亚型识别
MuSiC           - 组织细胞类型分解
EPIC            - 免疫细胞比例估计
MCPcounter      - 免疫细胞计数
ConsensusClusterPlus - 共识聚类
iClusterPlus    - 多组学数据整合聚类
```

## 🛠️ **开发与工具包**
```
devtools        - 包开发工具
usethis         - 自动化包设置
testthat        - 单元测试
pkgdown         - 包文档网站
roxygen2        - 文档生成
Rcpp            - R与C++集成
reticulate      - Python接口
```

## 📦 **基础与核心包**
```
Biobase         - 生物实验数据基础
BiocGenerics    - Bioconductor通用函数
BiocManager     - Bioconductor包管理
BiocParallel    - 并行计算
Matrix          - 稀疏矩阵处理
utils           - 工具函数
tools           - 包开发工具
stats           - 统计函数
graphics        - 基础图形
grDevices       - 图形设备
```
