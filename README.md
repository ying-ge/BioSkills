# FigureYa-package

600+个专业生物医学R包

这个列表几乎涵盖了现代生物医学研究的所有主要领域，从基础的数据处理到前沿的单细胞和空间组学技术。Bioconductor平台本身就有3000+个包，加上CRAN中的生物医学相关包，构成了极其丰富的生态系统。

这个生态系统仍在快速发展，新的包和功能不断涌现，特别是在单细胞、空间组学、机器学习和多组学整合领域。

## 🧬 **非编码RNA分析**
```
miRNAtap        - miRNA靶基因预测
multiMiR        - miRNA-mRNA互作数据库
mirBaseConverter - miRNA ID转换
miRNAtap        - 多数据库miRNA靶基因预测
lncRNAMap       - lncRNA功能注释
LncRNA2Target   - lncRNA靶基因数据库
LncRNADisease   - lncRNA疾病关联
circRNAprofiler - circRNA分析
circAtlas       - circRNA数据库
tRNA            - tRNA分析
```

## 🧪 **表观遗传学扩展**
```
ChIPseeker      - ChIP-seq峰注释
DiffBind        - ChIP-seq差异分析
csaw            - ChIP-seq峰检测
SICtools        - 表观遗传学工具
IlluminaHumanMethylationEPICanno.ilmn10.hg4 - EPIC芯片注释
EpiDISH         - 表观遗传去卷积
MethylMix       - 甲基化混合模型
```

## 🩸 **免疫组库分析**
```
immunarch       - 免疫组库分析核心包
VDJtools        - 免疫组库分析
Immunarch       - TCR/BCR分析
ChangeO         - 免疫组库分析
shazam          - 免疫组库统计分析
Alakazam        - 免疫组库可视化
```

## 🎯 **单细胞分析扩展包**
```
scater          - 单细胞质量控制
scran           - 单细胞标准化
scran           - 单细胞基因表达分析
scater          - 单细胞实验数据
scRNA-tools生态包
scater          - 单细胞数据可视化
scran           - 单细胞归一化
scater          - 单细胞特征选择
scran           - 单细胞差异表达
scater          - 单细胞降维
scran           - 单细胞聚类
scater          - 单细胞轨迹分析
scater          - 单细胞质量评估
scran           - 单细胞预处理
scater          - 单细胞实验设计
scran           - 单细胞技术对比
```

## 🎯 **深度学习在生物医学中的应用**
```
keras           - 深度学习框架
tensorflow      - 深度学习后端
torch           - PyTorch R接口
MXNet           - 深度学习框架
h2o             - 机器学习平台
caretEnsemble   - 集成学习
SuperLearner    - 超级学习器
mlr3            - 机器学习框架生态
tidymodels      - 整洁建模框架
```

## 🧪 **蛋白质结构和功能**
```
bio3d           - 生物分子3D结构
ProTraS         - 蛋白质翻译分析
protr           - 蛋白质特征计算
Rcintr          - 蛋白质相互作用中心
protViz         - 蛋白质组学可视化
Rdisop          - 蛋白质溶解度预测
```

## 🩺 **临床决策支持**
```
epic            - 内科病人预后预测
PredictABEL     - 疾病风险评估
riskRegression  - 风险回归分析
dynpred         - 动态预测模型
pec             - 预测误差曲线
```

## 🦠 **病毒学和传染病**
```
viral           - 病毒基因组分析
virusclust      - 病毒聚类分析
phyloscan       - 病毒系统发育扫描
OutbreakTools   - 疫情爆发分析
EpiEstim        - 疫情传播估计
EpiModel        - 传染病模型
```

## 🧬 **系统生物学扩展**
```
CellDesigner    - 细胞建模
synergyfinder   - 药物协同分析
SynergyFinder   - 药物协同效应
cellHTS2        - 高通量筛选分析
FlowSOM         - 流式数据自组织映射
Rcy3            - Cytoscape R接口
RCytoscape      - Cytoscape R接口
```

## 📊 **更多可视化包**
```
gganimate       - 动画图表
rayshader       - 3D可视化
rgl             - 3D图形
plotrix         - 特殊图表
ggforce         - 几何图形扩展
ggalt           - 替代坐标图
ggExtra         - 边际图
ggformula       - 公式接口
ggstance        - 位置调整
ggrepel         - 标签防重复
ggsci           - 科学配色
ggthemes        - 主题集合
ggtech          - 技术公司主题
```

## 🔬 **实验设计**
```
pwr             - 功效分析
pwr2            - 功效分析扩展
faraway         - 实验设计
agricolae       - 农业实验设计
DoE.base        - 实验设计基础
AlgDesign       - 实验设计算法
```

## 🧮 **生物统计学扩展**
```
multcomp        - 多重比较
lme4            - 线性混合模型
nlme            - 非线性混合模型
MCMCglmm        - 贝叶斯混合模型
brms            - 贝叶斯回归模型
rstanarm        - 贝叶斯回归
bayesplot       - 贝叶斯可视化
```

## 📈 **时间序列分析**
```
zoo             - 时间序列
xts             - 扩展时间序列
tseries         - 时间序列分析
forecast        - 时间序列预测
timetk          - 时间序列工具包
```

## 🎯 **临床数据管理**
```
clinicalutils   - 临床工具
mediana         - 临床试验设计
clinfun         - 临床函数
clinPK          - 临床药代动力学
clinReport      - 临床报告
```

## 🧪 **质谱数据扩展**
```
MSnbase         - 质谱数据处理
MSstatsTMT      - TMT标记定量
MSnID           - 质谱ID鉴定
MSqRob          - 质谱稳健统计
MZmine          - 质谱数据处理
OpenMS          - 开放质谱系统
```

## 🌐 **网络分析扩展**
```
igraphdata      - 图数据集
igraph         - 图分析(扩展包)
network         - 网络分析
sna            - 社会网络分析
tnet           - 时间网络分析
multiplex      - 多层网络
```

## 🧬 **基因组浏览器**
```
Gviz           - 基因组可视化
trackViewer    - 轨道查看器
ggbio          - ggplot2基因组
GenomicRanges  - 基因组区间(扩展)
rtracklayer    - 基因组轨道
```

## 🎯 **多组学整合**
```
MOFA           - 多组学因子分析
mixOmics       - 多组学分析
iClusterPlus   - 整合聚类
SNFtool        - 相似性网络融合
MultiAssayExperiment - 多组学实验
DIABLO         - 多组学数据整合
```

## 🔬 **更多机器学习包**
```
xgboost        - 梯度提升
lightgbm       - 轻量级梯度提升
catboost       - 分类提升
randomForestSRC - 随机森林生存
ranger         - 快速随机森林
mlr3learners   - mlr3学习器
tidymodels     - 整洁建模
parsnip        - 模型接口
workflows      - 工作流
tune           - 超参数调优
```

## 📊 **大数据处理**
```
data.table     - 大数据处理(扩展)
disk.frame     - 磁盘数据框
bigmemory      - 大内存对象
ff             - 文件因子
bigstatsr      - 大统计
biglasso       - 大数据lasso
```

## 🎯 **空间组学扩展**
```
spacexr        - 空间表达数据
spatstat       - 空间统计
sf             - 简单特征
sp             - 空间数据
raster         - 栅格数据
terra          - 现代栅格数据
```

## 🧬 **新兴技术包**
```
squidpy        - 空间转录组Python接口(通过reticulate)
scVI           - 单细胞变分推断
scANVI         - 单细胞注释变分推断
totalVI        - 总RNA单细胞分析
```

## 📱 **报告和文档**
```
rmarkdown      - 动态文档
knitr          - 动态报告
bookdown       - 书籍编写
blogdown       - 博客系统
shiny          - 交互式应用
flexdashboard  - 灵活仪表板
```

## 🎯 **生物信息学工作流**
```
drake          - 数据管道
targets        - 目标管道
workflowr      - 工作流管理
renv           - 环境管理
packrat        - 包管理
```
