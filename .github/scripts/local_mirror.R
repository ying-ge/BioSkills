# 创建本地CRAN镜像风格的仓库
create_local_repo <- function(backup_dir = "CRAN_Bioc_backup", repo_dir = "local_repo") {
  if (!requireNamespace("tools", quietly = TRUE)) {
    stop("需要tools包")
  }
  
  # 创建仓库目录结构
  src_contrib <- file.path(repo_dir, "src", "contrib")
  dir.create(src_contrib, recursive = TRUE)
  
  # 复制所有tar.gz文件
  tar_files <- list.files(backup_dir, pattern = "\\.tar\\.gz$", 
                         recursive = TRUE, full.names = TRUE)
  
  for (tar_file in tar_files) {
    file.copy(tar_file, src_contrib)
    cat("复制:", basename(tar_file), "\n")
  }
  
  # 生成PACKAGES文件
  setwd(src_contrib)
  tools::write_PACKAGES(".", type = "source")
  setwd("../../..")
  
  cat("✅ 本地仓库创建完成:", normalizePath(repo_dir), "\n")
  cat("使用方法:\n")
  cat('options(repos = c(LOCAL = "file://', normalizePath(repo_dir), '", CRAN = "https://cloud.r-project.org/"))\n')
}
 
# 使用本地仓库
use_local_repo <- function(repo_dir = "local_repo") {
  local_repo_url <- paste0("file://", normalizePath(repo_dir))
  options(repos = c(LOCAL = local_repo_url, getOption("repos")))
  cat("✅ 本地仓库已添加到repos选项\n")
  cat("现在可以使用 install.packages() 从本地仓库安装\n")
}
