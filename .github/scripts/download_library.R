# 设置下载目录
download_dir <- .
if (!dir.exists(download_dir)) {
  dir.create(download_dir)
}

# 包列表
packages <- c(
  "AUCell", "AnnotationDbi", "AnnotationHub", "BART", "BSgenome", 
  "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", 
  "Biobase", "BiocParallel", "BlandAltmanLeh", "Boruta", "CMScaller", 
  "Cairo", "CancerSubtypes", "CellChat", "ChAMP", "ChAMPdata", 
  "ChIPseeker", "ClassDiscovery", "ComparisonSurv", "ComplexHeatmap", 
  "ConsensusClusterPlus", "CoxBoost", "DESeq2", "DNAcopy", "DOSE", 
  "DT", "DealGPL570", "EPIC", "EnvStats", "GEOquery", "GO.db", 
  "GOSemSim", "GOplot", "GSA", "GSEABase", "GSVA", "GenomicFeatures", 
  "GenomicRanges", "GetoptLong", "Gviz", "Hmisc", "IOBR", "ISOpureR", 
  "IlluminaHumanMethylation450kanno.ilmn12.hg19", "ImmuLncRNA", 
  "KEGGREST", "MCPcounter", "MatchIt", "Matrix", "MethComp", "MuSiC", 
  "MutationalPatterns", "NMF", "OmicCircos", "PMAPscore", "PharmacoGx", 
  "RCircos", "RColorBrewer", "RIdeogram", "RTN", "ResourceSelection", 
  "Rmisc", "S4Vectors", "SNFtool", "SNPlocs.Hsapiens.dbSNP142.GRCh37", 
  "SPIA", "ScBulkCCCI", "Seurat", "SeuratData", "SeuratDisk", 
  "SimDesign", "SingleCellExperiment", "SingleR", "SummarizedExperiment", 
  "SuperExactTest", "Sushi", "TCGAbiolinks", 
  "TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
  "VennDiagram", "WGCNA", "affy", "affyPLM", "ape", "aplot", "bamp", 
  "biomaRt", "caTools", "car", "caret", "cgdsr", "cicero", "circlize", 
  "classpredict", "cluster", "clusterProfiler", "coin", "colorRamps", 
  "compareC", "condsurv", "corrplot", "cowplot", "data.table", 
  "dendsort", "devtools", "directlabels", "doParallel", "dorothea", 
  "dplyr", "e1071", "edgeR", "enrichplot", "entropy", "epitools", 
  "estimate", "export", "extrafont", "factoextra", "fgsea", "fishplot", 
  "fmsb", "forcats", "foreign", "forestplot", "future", "future.apply", 
  "gbm", "gcrma", "genefu", "ggVennDiagram", "ggalluvial", "gganatogram", 
  "ggbreak", "ggcor", "ggcorrplot", "ggnewscale", "ggpattern", "ggplot2", 
  "ggplotify", "ggpp", "ggprism", "ggpubr", "ggraph", "ggrepel", 
  "ggridges", "ggsci", "ggsignif", "ggtext", "ggthemes", "ggtree", 
  "givitiR", "glmGamPoi", "glmnet", "glue", "gplots", "grDevices", 
  "graphics", "grid", "gridBase", "gridExtra", "gtable", "gtools", 
  "hgu133plus2.db", "htmlwidgets", "iClusterPlus", "ifnb.SeuratData", 
  "igraph", "immunedeconv", "impute", "infotheo", "janitor", "jsonlite", 
  "kableExtra", "ktplots", "lattice", "limma", "lme4", "lolR", "mRMRe", 
  "maftools", "magick", "magrittr", "matlab", "matrixStats", "mclust", 
  "meta", "miRBaseConverter", "miRBaseVersions.db", "minfi", "miscTools", 
  "mixtools", "mlbench", "monocle", "monocle3", "motifbreakR", "msa", 
  "msigdbr", "muhaz", "networkD3", "nlme", "nomogramEx", "officer", 
  "oncoPredict", "oompaBase", "openxlsx", "org.Hs.eg.db", "org.Mm.eg.db", 
  "pROC", "pamr", "parallel", "patchwork", "pathifier", "pbapply", 
  "pec", "phangorn", "pheatmap", "plsRcox", "plyr", "preprocessCore", 
  "presto", "pubmed.mineR", "purrr", "randomForest", "randomForestSRC", 
  "randomSurvivalForest", "randomcoloR", "ranger", "readr", "readxl", 
  "regplot", "rentrez", "reshape2", "rgl", "rhdf5", "ridge", "rlang", 
  "rms", "rstatix", "rsvd", "rsvg", "rtracklayer", "scMetabolism", 
  "scales", "seqinr", "shape", "showtext", "sigFeature", "slingshot", 
  "smatr", "smoothHR", "snow", "spatstat.explore", "squash", "stringi", 
  "stringr", "strucchange", "superpc", "survMisc", "survRM2", "survcomp", 
  "survival", "survivalROC", "survivalsvm", "survminer", "sva", "tibble", 
  "tidyHeatmap", "tidyestimate", "tidygraph", "tidyjson", "tidyr", 
  "tidyverse", "timeROC", "tools", "trackViewer", "tradeSeq", "usedist", 
  "utils", "vegan", "viper", "viridis", "wakefield", "wateRmelon", 
  "xgboost", "xlsx"
)

# 初始化结果记录
results <- data.frame(
  package = packages,
  status = "未尝试",
  source = NA,
  file_path = NA,
  error_message = NA,
  stringsAsFactors = FALSE
)

# 设置下载选项
options(timeout = 300) # 设置超时时间为5分钟

# 检查并安装必要的包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# 下载函数
download_package <- function(pkg) {
  cat("正在下载:", pkg, "\n")
  
  # 尝试从CRAN下载
  try_cran <- tryCatch({
    cran_url <- paste0("https://cran.r-project.org/src/contrib/", pkg, "_*.tar.gz")
    available_pkgs <- available.packages()
    if (pkg %in% rownames(available_pkgs)) {
      download.packages(pkg, destdir = download_dir, type = "source")
      return(list(status = "成功", source = "CRAN"))
    } else {
      return(list(status = "失败", source = "CRAN"))
    }
  }, error = function(e) {
    return(list(status = "失败", source = "CRAN", error = e$message))
  })
  
  if (try_cran$status == "成功") {
    return(try_cran)
  }
  
  # 尝试从Bioconductor下载
  try_bioc <- tryCatch({
    bioc_url <- tryCatch({
      BiocManager::available(pkg)
      download.packages(pkg, destdir = download_dir, type = "source", repos = BiocManager::repositories())
      return(list(status = "成功", source = "Bioconductor"))
    }, error = function(e) {
      return(list(status = "失败", source = "Bioconductor", error = e$message))
    })
  }, error = function(e) {
    return(list(status = "失败", source = "Bioconductor", error = e$message))
  })
  
  if (try_bioc$status == "成功") {
    return(try_bioc)
  }
  
  # 尝试从GitHub下载（需要知道仓库名）
  try_github <- tryCatch({
    # 这里需要知道GitHub仓库名，我们可以尝试一些常见的模式
    # 注意：这只是一个示例，实际使用时需要更精确的映射
    if (pkg == "ScBulkCCCI") {
      devtools::install_github("ChenWeiRong/ScBulkCCCI", build = FALSE)
      return(list(status = "成功", source = "GitHub"))
    } else if (pkg == "ktplots") {
      devtools::install_github("zktuong/ktplots", build = FALSE)
      return(list(status = "成功", source = "GitHub"))
    } else {
      return(list(status = "失败", source = "GitHub", error = "未知的GitHub仓库"))
    }
  }, error = function(e) {
    return(list(status = "失败", source = "GitHub", error = e$message))
  })
  
  if (try_github$status == "成功") {
    return(try_github)
  }
  
  return(list(status = "失败", source = "所有源", error = "所有下载方法都失败"))
}

# 逐个下载包
for (i in 1:length(packages)) {
  pkg <- packages[i]
  cat("处理包 (", i, "/", length(packages), "):", pkg, "\n")
  
  result <- download_package(pkg)
  
  results$status[i] <- result$status
  results$source[i] <- result$source
  if (!is.null(result$error)) {
    results$error_message[i] <- result$error
  }
  
  # 添加延迟，避免请求过于频繁
  Sys.sleep(1)
}

# 保存结果
write.csv(results, file.path(download_dir, "download_results.csv"), row.names = FALSE)

# 输出摘要
success_count <- sum(results$status == "成功")
failed_count <- sum(results$status == "失败")

cat("\n=== 下载结果摘要 ===\n")
cat("成功下载:", success_count, "个包\n")
cat("下载失败:", failed_count, "个包\n")

cat("\n下载来源分布:\n")
print(table(results$source[results$status == "成功"]))

cat("\n失败的包:\n")
failed_packages <- results$package[results$status == "失败"]
if (length(failed_packages) > 0) {
  print(failed_packages)
} else {
  cat("无\n")
}

# 显示下载目录中的文件
cat("\n下载的文件列表:\n")
downloaded_files <- list.files(download_dir, pattern = "\\.tar\\.gz$")
print(downloaded_files)
