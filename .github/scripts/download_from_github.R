# 安装必要的工具包
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("packrat", quietly = TRUE)) {
  install.packages("packrat")
}
if (!requireNamespace("miniCRAN", quietly = TRUE)) {
  install.packages("miniCRAN")
}

library(devtools)
library(remotes)

# GitHub仓库列表
github_repos <- c(
  "tpq/kpmt",
  "cBioPortal/cgdsr",
  "satijalab/seurat-data",
  "jespermaag/gganatogram",
  "chrisamiller/fishplot",
  "BITS-VIB/venn-tools",
  "FredHutch/Galeano-Nino-Bullman-Intratumoral-Microbiota_2022",
  "Gibbsdavidl/Immune-Subtype-Clustering",
  "Gibbsdavidl/ImmuneSubtypeClassifier",
  "GuangchuangYu/gglayer",
  "IOBR/IOBR",
  "IndrajeetPatil/ggstatsplot",
  "KrishnaswamyLab/MAGIC",
  "NKI-CCB/DISCOVER",
  "Simon-Coetzee/motifBreakR",
  "Teichlab/cellphonedb",
  "YuLab-SMU/clusterProfiler-book",
  "aertslab/pySCENIC",
  "califano-lab/ARACNe-AP",
  "chrisamiller/fishplot",
  "cole-trapnell-lab/monocle-release",
  "cran/colorfulVennPlot",
  "dylkot/cNMF",
  "ebecht/MCPcounter",
  "fawda123/ggord",
  "icbi-lab/Immunophenogram",
  "inutano/pfastq-dump",
  "jespermaag/gganatogram",
  "johncolby/SVM-RFE",
  "linnarsson-lab/ipynb-lamanno2016",
  "lmcinnes/umap",
  "matanhofree/crc-immune-hubs",
  "matt-black/dcapy",
  "mg14/mg14",
  "nicholasehamilton/ggtern",
  "paulgeeleher/pRRophetic",
  "paulgeeleher/pRRophetic2",
  "peterawe/CMScaller",
  "raerose01/deconstructSigs",
  "rfordatascience/tidytuesday",
  "richie2019/MBpanel",
  "ropensci/rentrez",
  "rsankowski/sankowski-et-al-microglia",
  "spren9er/tidytuesday",
  "taoliu/MACS",
  "wanyewang1/Index_model",
  "woobe/rPlotter",
  "xlucpu/MOVICS",
  "zzwch/crosslink"
)

# 去除重复项
github_repos <- unique(github_repos)

# 创建下载目录
download_dir <- "r_packages_with_deps"
if (!dir.exists(download_dir)) {
  dir.create(download_dir)
}

# 函数：获取包的依赖关系
get_package_dependencies <- function(repo) {
  cat("分析依赖关系:", repo, "\n")
  
  tryCatch({
    # 获取包的描述信息
    desc_url <- paste0("https://raw.githubusercontent.com/", repo, "/master/DESCRIPTION")
    temp_file <- tempfile()
    download.file(desc_url, temp_file, quiet = TRUE)
    
    if (file.exists(temp_file)) {
      desc <- read.dcf(temp_file)
      deps <- c()
      
      # 获取Imports, Depends, LinkingTo字段
      dependency_fields <- c("Imports", "Depends", "LinkingTo")
      for (field in dependency_fields) {
        if (field %in% colnames(desc)) {
          field_deps <- desc[, field]
          # 解析依赖关系
          if (!is.na(field_deps)) {
            # 分割并清理依赖项
            pkgs <- strsplit(field_deps, ",")[[1]]
            pkgs <- gsub("\\([^)]*\\)", "", pkgs)  # 移除版本限制
            pkgs <- gsub("\\n", "", pkgs)  # 移除换行符
            pkgs <- trimws(pkgs)  # 移除空格
            pkgs <- pkgs[pkgs != "R"]  # 移除R本身
            deps <- c(deps, pkgs)
          }
        }
      }
      
      file.remove(temp_file)
      return(unique(deps))
    }
    return(character(0))
  }, error = function(e) {
    cat("无法获取", repo, "的依赖关系:", e$message, "\n")
    return(character(0))
  })
}

# 函数：下载CRAN包及其依赖
download_cran_package_with_deps <- function(pkg_name, download_dir) {
  cat("下载CRAN包及依赖:", pkg_name, "\n")
  
  tryCatch({
    # 获取包的所有依赖
    deps <- tools::package_dependencies(pkg_name, recursive = TRUE, 
                                       which = c("Depends", "Imports", "LinkingTo"))
    all_pkgs <- unique(c(pkg_name, unlist(deps)))
    
    # 下载所有包
    downloaded <- download.packages(all_pkgs, destdir = download_dir, type = "source")
    
    return(list(status = "成功", packages = downloaded[,1], files = downloaded[,2]))
  }, error = function(e) {
    return(list(status = "失败", error = e$message))
  })
}

# 函数：安装GitHub包并处理依赖
install_github_with_deps <- function(repo, download_dir) {
  pkg_name <- strsplit(repo, "/")[[1]][2]
  
  cat("\n=== 处理包:", repo, "===\n")
  
  # 1. 首先获取依赖关系
  dependencies <- get_package_dependencies(repo)
  cat("发现的依赖:", if (length(dependencies) > 0) paste(dependencies, collapse = ", ") else "无", "\n")
  
  # 2. 下载所有CRAN依赖
  if (length(dependencies) > 0) {
    for (dep in dependencies) {
      cat("下载依赖:", dep, "\n")
      download_cran_package_with_deps(dep, download_dir)
    }
  }
  
  # 3. 下载GitHub包源代码
  cat("下载GitHub包:", repo, "\n")
  tryCatch({
    # 尝试main分支
    download.file(
      paste0("https://github.com/", repo, "/archive/main.tar.gz"),
      destfile = file.path(download_dir, paste0(pkg_name, "_github.tar.gz")),
      quiet = FALSE
    )
    github_status <- "成功"
  }, error = function(e) {
    # 尝试master分支
    tryCatch({
      download.file(
        paste0("https://github.com/", repo, "/archive/master.tar.gz"),
        destfile = file.path(download_dir, paste0(pkg_name, "_github.tar.gz")),
        quiet = FALSE
      )
      github_status <- "成功"
    }, error = function(e2) {
      github_status <- "失败"
    })
  })
  
  return(list(
    package = pkg_name,
    repository = repo,
    dependencies = dependencies,
    github_download_status = github_status
  ))
}

# 批量下载所有包及其依赖
all_results <- list()

for (i in 1:length(github_repos)) {
  repo <- github_repos[i]
  result <- install_github_with_deps(repo, download_dir)
  all_results[[i]] <- result
  
  # 添加延迟
  Sys.sleep(2)
}

# 汇总结果
summary_df <- do.call(rbind, lapply(all_results, function(x) {
  data.frame(
    package = x$package,
    repository = x$repository,
    dependencies_count = length(x$dependencies),
    dependencies = paste(x$dependencies, collapse = ", "),
    github_download_status = x$github_download_status,
    stringsAsFactors = FALSE
  )
}))

# 保存结果
write.csv(summary_df, file.path(download_dir, "dependency_analysis_results.csv"), row.names = FALSE)

# 下载常见的CRAN依赖包（预先下载）
common_cran_deps <- c(
  "ggplot2", "dplyr", "tidyr", "purrr", "stringr", "readr",
  "data.table", "magrittr", "tibble", "httr", "jsonlite",
  "rlang", "glue", "withr", "R6", "methods", "utils", "stats",
  "graphics", "grDevices", "grid", "parallel", "future",
  "doParallel", "foreach", "iterators", "BiocManager", "remotes",
  "devtools", "testthat", "knitr", "rmarkdown", "roxygen2"
)

cat("\n=== 下载常见的CRAN依赖包 ===\n")
for (pkg in common_cran_deps) {
  if (!pkg %in% rownames(installed.packages())) {
    cat("下载:", pkg, "\n")
    download.packages(pkg, destdir = download_dir, type = "source")
  }
}

# 显示下载的文件
cat("\n=== 下载的文件列表 ===\n")
downloaded_files <- list.files(download_dir, pattern = "\\.tar\\.gz$")
cat("总共下载了", length(downloaded_files), "个文件\n")
print(downloaded_files)

# 生成依赖关系报告
cat("\n=== 依赖关系摘要 ===\n")
cat("总共处理了", nrow(summary_df), "个GitHub包\n")
cat("平均每个包有", mean(summary_df$dependencies_count), "个依赖\n")

# 显示依赖最多的包
if (nrow(summary_df) > 0) {
  most_deps <- summary_df[order(-summary_df$dependencies_count), ][1:5, ]
  cat("\n依赖最多的前5个包:\n")
  print(most_deps[, c("package", "dependencies_count")])
}

# 保存会话信息
sink(file.path(download_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== 下载完成 ===\n")
cat("所有包和依赖已下载到:", download_dir, "\n")
cat("详细结果保存在: dependency_analysis_results.csv\n")
