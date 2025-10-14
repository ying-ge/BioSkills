# 从GitHub备份安装R包
install_from_github_backup <- function(repo_name, backup_dir = "R_packages_backup") {
  # 转换仓库名为目录名
  repo_dir_name <- gsub("/", "_", repo_name)
  repo_path <- file.path(backup_dir, repo_dir_name)
  
  if (!dir.exists(repo_path)) {
    stop(paste("备份目录不存在:", repo_path))
  }
  
  # 查找zip文件
  zip_files <- list.files(repo_path, pattern = "\\.zip$", full.names = TRUE)
  if (length(zip_files) == 0) {
    stop(paste("在", repo_path, "中未找到zip文件"))
  }
  
  # 使用第一个zip文件
  zip_file <- zip_files[1]
  cat("正在安装:", repo_name, "从", basename(zip_file), "\n")
  
  # 解压到临时目录
  temp_dir <- tempdir()
  extract_dir <- file.path(temp_dir, paste0("extract_", Sys.time() %>% as.numeric()))
  dir.create(extract_dir)
  
  unzip(zip_file, exdir = extract_dir)
  
  # 找到解压后的包目录
  extracted_dirs <- list.dirs(extract_dir, recursive = FALSE)
  if (length(extracted_dirs) == 0) {
    stop("解压失败或未找到包目录")
  }
  
  pkg_dir <- extracted_dirs[1]
  
  # 检查是否是有效的R包
  desc_file <- file.path(pkg_dir, "DESCRIPTION")
  if (!file.exists(desc_file)) {
    stop("这不是一个有效的R包目录")
  }
  
  # 安装包
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  
  devtools::install(pkg_dir, dependencies = TRUE)
  
  # 清理临时文件
  unlink(extract_dir, recursive = TRUE)
  
  cat("✅", repo_name, "安装完成\n")
}
 
# 使用示例
install_from_github_backup("tpq/kpmt")
install_from_github_backup("ggplot2/ggplot2")
