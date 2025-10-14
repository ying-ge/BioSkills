# 查看GitHub备份信息
list_github_backups <- function(backup_dir = "R_packages_backup") {
  if (!dir.exists(backup_dir)) {
    cat("备份目录不存在\n")
    return(invisible(NULL))
  }
  
  # 读取日志文件
  log_files <- list.files(backup_dir, pattern = "backup_log.*\\.csv$", full.names = TRUE)
  if (length(log_files) > 0) {
    log_data <- read.csv(log_files[1])
    cat("GitHub备份统计:\n")
    print(table(log_data$status))
    
    successful <- log_data[log_data$status == "success", ]
    if (nrow(successful) > 0) {
      cat("\n成功备份的仓库:\n")
      print(successful[, c("repo", "file_size_mb")])
    }
  }
  
  # 列出所有备份目录
  backup_dirs <- list.dirs(backup_dir, recursive = FALSE)
  cat("\n备份目录数量:", length(backup_dirs), "\n")
}
 
# 查看CRAN/Bioconductor备份信息
list_cran_backups <- function(backup_dir = "CRAN_Bioc_backup") {
  if (!dir.exists(backup_dir)) {
    cat("备份目录不存在\n")
    return(invisible(NULL))
  }
  
  # 读取日志文件
  log_files <- list.files(backup_dir, pattern = "backup_log.*\\.csv$", full.names = TRUE)
  if (length(log_files) > 0) {
    log_data <- read.csv(log_files[1])
    cat("CRAN/Bioconductor备份统计:\n")
    print(table(log_data$status))
    print(table(log_data$source))
    
    successful <- log_data[log_data$status == "success", ]
    if (nrow(successful) > 0) {
      cat("\n成功备份的包:\n")
      print(successful[, c("package", "version", "source", "file_size_mb")])
    }
  }
  
  # 列出所有备份目录
  backup_dirs <- list.dirs(backup_dir, recursive = FALSE)
  cat("\n备份包数量:", length(backup_dirs), "\n")
}
