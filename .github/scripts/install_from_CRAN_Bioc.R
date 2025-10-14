# 从CRAN/Bioconductor备份安装R包
install_from_cran_backup <- function(package_name, backup_dir = "CRAN_Bioc_backup") {
  pkg_path <- file.path(backup_dir, package_name)
  
  if (!dir.exists(pkg_path)) {
    stop(paste("备份目录不存在:", pkg_path))
  }
  
  # 查找tar.gz文件
  tar_files <- list.files(pkg_path, pattern = "\\.tar\\.gz$", full.names = TRUE)
  if (length(tar_files) == 0) {
    stop(paste("在", pkg_path, "中未找到tar.gz文件"))
  }
  
  # 使用第一个tar.gz文件
  tar_file <- tar_files[1]
  cat("正在安装:", package_name, "从", basename(tar_file), "\n")
  
  # 直接从tar.gz文件安装
  install.packages(tar_file, repos = NULL, type = "source")
  
  cat("✅", package_name, "安装完成\n")
}
 
# 使用示例
install_from_cran_backup("ggplot2")
install_from_cran_backup("limma")
