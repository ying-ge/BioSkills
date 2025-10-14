# scripts/download_packages.R
#!/usr/bin/env Rscript
# 下载R包到本地仓库
 
suppressPackageStartupMessages({
  library(jsonlite)
  library(utils)
})
 
PACKAGE_INPUT <- "metadata/figureya_packages.json"
PACKAGE_DIR <- "packages"
LOG_DIR <- "metadata/download_logs"
 
# 创建目录结构
setup_directories <- function() {
  dir.create(PACKAGE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PACKAGE_DIR, "CRAN"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PACKAGE_DIR, "Bioconductor"), showWarnings = FALSE, recursive = TRUE)
  dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
}
 
# 下载CRAN包
download_cran_package <- function(package_name, version, r_version) {
  r_major_minor <- sub("\\.[0-9]+$", "", r_version)
  pkg_dir <- file.path(PACKAGE_DIR, "CRAN", r_major_minor)
  dir.create(pkg_dir, showWarnings = FALSE, recursive = TRUE)
  
  pkg_file <- paste0(package_name, "_", version, ".tar.gz")
  pkg_path <- file.path(pkg_dir, pkg_file)
  
  if(file.exists(pkg_path)) {
    return(list(success = TRUE, message = "Already exists", path = pkg_path))
  }
  
  # 尝试多个下载源
  urls <- c(
    paste0("https://cran.r-project.org/src/contrib/Archive/", package_name, "/", pkg_file),
    paste0("https://cran.r-project.org/src/contrib/", pkg_file),
    paste0("https://cloud.r-project.org/src/contrib/Archive/", package_name, "/", pkg_file),
    paste0("https://cloud.r-project.org/src/contrib/", pkg_file)
  )
  
  for(url in urls) {
    result <- tryCatch({
      download.file(url, pkg_path, mode = "wb", quiet = TRUE)
      if(file.exists(pkg_path) && file.size(pkg_path) > 0) {
        return(list(success = TRUE, message = paste("Downloaded from", url), path = pkg_path))
      }
      FALSE
    }, error = function(e) FALSE)
    
    if(isTRUE(result$success)) break
  }
  
  return(list(success = FALSE, message = "All download attempts failed", path = NA))
}
 
# 主下载函数
main <- function() {
  setup_directories()
  
  if(!file.exists(PACKAGE_INPUT)) {
    stop("Package input file not found: ", PACKAGE_INPUT)
  }
  
  cat("读取包信息...\n")
  data <- fromJSON(PACKAGE_INPUT)
  
  download_log <- list(
    start_time = Sys.time(),
    total_packages = length(data$packages),
    results = list(),
    summary = list(success = 0, failed = 0, skipped = 0)
  )
  
  cat("开始下载", length(data$packages), "个包...\n")
  
  for(pkg_name in names(data$packages)) {
    cat("处理包:", pkg_name, "\n")
    
    pkg_versions <- data$packages[[pkg_name]]
    
    for(r_version in names(pkg_versions)) {
      pkg_info <- pkg_versions[[r_version]]
      
      if(pkg_info$source == "CRAN") {
        result <- download_cran_package(pkg_name, pkg_info$version, r_version)
        
        download_log$results[[paste(pkg_name, r_version, sep = "_")]] <- list(
          package = pkg_name,
          version = pkg_info$version,
          r_version = r_version,
          source = pkg_info$source,
          success = result$success,
          message = result$message,
          path = result$path,
          timestamp = Sys.time()
        )
        
        if(result$success) {
          download_log$summary$success <- download_log$summary$success + 1
          cat("  ✓", pkg_name, pkg_info$version, "(R", r_version, ")\n")
        } else {
          download_log$summary$failed <- download_log$summary$failed + 1
          cat("  ✗", pkg_name, pkg_info$version, "(R", r_version, "):", result$message, "\n")
        }
      } else {
        download_log$summary$skipped <- download_log$summary$skipped + 1
        cat("  - 跳过非CRAN包:", pkg_name, "\n")
      }
    }
  }
  
  download_log$end_time <- Sys.time()
  download_log$duration <- difftime(download_log$end_time, download_log$start_time, units = "mins")
  
  # 保存下载日志
  log_file <- file.path(LOG_DIR, paste0("download_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
  writeLines(toJSON(download_log, pretty = TRUE, auto_unbox = TRUE), log_file)
  
  cat("\n下载完成!\n")
  cat("成功:", download_log$summary$success, "\n")
  cat("失败:", download_log$summary$failed, "\n")
  cat("跳过:", download_log$summary$skipped, "\n")
  cat("耗时:", round(as.numeric(download_log$duration), 2), "分钟\n")
  cat("日志文件:", log_file, "\n")
}
 
if(!interactive()) {
  main()
}
