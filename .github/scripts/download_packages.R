# scripts/download_packages_with_size_check.R
#!/usr/bin/env Rscript
# 下载R包到本地仓库，带文件大小检查
 
suppressPackageStartupMessages({
  library(jsonlite)
  library(utils)
})
 
PACKAGE_INPUT <- "metadata/figureya_packages.json"
PACKAGE_DIR <- "packages"
LOG_DIR <- "metadata/download_logs"
MAX_FILE_SIZE_MB <- 100  # GitHub限制
 
# 创建目录结构
setup_directories <- function() {
  dir.create(PACKAGE_DIR, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PACKAGE_DIR, "CRAN"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(PACKAGE_DIR, "Bioconductor"), showWarnings = FALSE, recursive = TRUE)
  dir.create(LOG_DIR, showWarnings = FALSE, recursive = TRUE)
}
 
# 检查文件大小
check_file_size <- function(file_path, max_size_mb = MAX_FILE_SIZE_MB) {
  if(!file.exists(file_path)) return(list(valid = FALSE, size_mb = 0))
  
  size_bytes <- file.size(file_path)
  size_mb <- size_bytes / (1024 * 1024)
  
  return(list(
    valid = size_mb <= max_size_mb,
    size_mb = round(size_mb, 2),
    size_bytes = size_bytes
  ))
}
 
# 下载CRAN包
download_cran_package <- function(package_name, version, r_version) {
  r_info <- R.version
  r_major_minor <- paste(r_info$major, strsplit(r_info$minor, "\\.")[[1]][1], sep = ".")
  
  pkg_dir <- file.path(PACKAGE_DIR, "CRAN", r_major_minor)
  dir.create(pkg_dir, showWarnings = FALSE, recursive = TRUE)
  
  pkg_file <- paste0(package_name, "_", version, ".tar.gz")
  pkg_path <- file.path(pkg_dir, pkg_file)
  
  if(file.exists(pkg_path)) {
    size_check <- check_file_size(pkg_path)
    return(list(
      success = TRUE, 
      message = "Already exists", 
      path = pkg_path,
      size_mb = size_check$size_mb,
      size_valid = size_check$valid
    ))
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
        size_check <- check_file_size(pkg_path)
        
        if(!size_check$valid) {
          # 文件太大，移动到特殊目录
          large_file_path <- paste0(pkg_path, ".large")
          file.rename(pkg_path, large_file_path)
          
          return(list(
            success = TRUE, 
            message = paste("Downloaded but too large (", size_check$size_mb, "MB) - moved to .large"),
            path = large_file_path,
            size_mb = size_check$size_mb,
            size_valid = FALSE,
            url = url
          ))
        }
        
        return(list(
          success = TRUE, 
          message = paste("Downloaded from", url),
          path = pkg_path,
          size_mb = size_check$size_mb,
          size_valid = TRUE,
          url = url
        ))
      }
      FALSE
    }, error = function(e) FALSE)
    
    if(isTRUE(result$success)) break
  }
  
  return(list(
    success = FALSE, 
    message = "All download attempts failed", 
    path = NA,
    size_mb = 0,
    size_valid = FALSE
  ))
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
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    max_file_size_mb = MAX_FILE_SIZE_MB,
    results = list(),
    summary = list(
      success = 0, 
      failed = 0, 
      skipped = 0, 
      too_large = 0,
      total_size_mb = 0
    )
  )
  
  cat("开始下载", length(data$packages), "个包...\n")
  cat("当前R版本:", download_log$r_version, "\n")
  cat("文件大小限制:", MAX_FILE_SIZE_MB, "MB\n\n")
  
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
          size_mb = result$size_mb,
          size_valid = result$size_valid,
          timestamp = Sys.time()
        )
        
        if(result$success) {
          if(result$size_valid) {
            download_log$summary$success <- download_log$summary$success + 1
            download_log$summary$total_size_mb <- download_log$summary$total_size_mb + result$size_mb
            cat("  ✅", pkg_name, pkg_info$version, sprintf("(%.1fMB)", result$size_mb), "\n")
          } else {
            download_log$summary$too_large <- download_log$summary$too_large + 1
            cat("  ⚠️ ", pkg_name, pkg_info$version, sprintf("(%.1fMB - TOO LARGE)", result$size_mb), "\n")
          }
        } else {
          download_log$summary$failed <- download_log$summary$failed + 1
          cat("  ❌", pkg_name, pkg_info$version, ":", result$message, "\n")
        }
      } else {
        download_log$summary$skipped <- download_log$summary$skipped + 1
        cat("  ⏭️  跳过非CRAN包:", pkg_name, "\n")
      }
    }
  }
  
  download_log$end_time <- Sys.time()
  download_log$duration <- difftime(download_log$end_time, download_log$start_time, units = "mins")
  download_log$summary$total_size_mb <- round(download_log$summary$total_size_mb, 2)
  
  # 保存下载日志
  log_file <- file.path(LOG_DIR, paste0("download_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".json"))
  writeLines(toJSON(download_log, pretty = TRUE, auto_unbox = TRUE), log_file)
  
  cat("\n📊 下载完成!\n")
  cat("✅ 成功:", download_log$summary$success, "\n")
  cat("❌ 失败:", download_log$summary$failed, "\n")
  cat("⏭️  跳过:", download_log$summary$skipped, "\n")
  cat("⚠️  超大文件:", download_log$summary$too_large, "\n")
  cat("📦 总大小:", download_log$summary$total_size_mb, "MB\n")
  cat("⏱️  耗时:", round(as.numeric(download_log$duration), 2), "分钟\n")
  cat("📄 日志文件:", log_file, "\n")
  
  # 如果有超大文件，给出建议
  if(download_log$summary$too_large > 0) {
    cat("\n💡 建议:\n")
    cat("   - 超大文件已标记为.large，不会提交到git\n")
    cat("   - 考虑将这些文件上传到GitHub Releases\n")
    cat("   - 或使用Git LFS进行大文件管理\n")
  }
}
 
if(!interactive()) {
  main()
}
