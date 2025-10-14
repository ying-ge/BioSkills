# 批量安装GitHub备份包
batch_install_github_backup <- function(repo_list, backup_dir = "R_packages_backup") {
  success_count <- 0
  failed_list <- character()
  
  for (repo in repo_list) {
    cat("\n正在安装:", repo, "\n")
    tryCatch({
      install_from_github_backup(repo, backup_dir)
      success_count <- success_count + 1
    }, error = function(e) {
      cat("❌ 安装失败:", repo, "-", e$message, "\n")
      failed_list <<- c(failed_list, repo)
    })
  }
  
  cat("\n=== 安装完成 ===\n")
  cat("成功:", success_count, "\n")
  cat("失败:", length(failed_list), "\n")
  if (length(failed_list) > 0) {
    cat("失败列表:\n", paste(failed_list, collapse = "\n"), "\n")
  }
}
 
# 批量安装CRAN/Bioconductor备份包
batch_install_cran_backup <- function(package_list, backup_dir = "CRAN_Bioc_backup") {
  success_count <- 0
  failed_list <- character()
  
  for (pkg in package_list) {
    cat("\n正在安装:", pkg, "\n")
    tryCatch({
      install_from_cran_backup(pkg, backup_dir)
      success_count <- success_count + 1
    }, error = function(e) {
      cat("❌ 安装失败:", pkg, "-", e$message, "\n")
      failed_list <<- c(failed_list, pkg)
    })
  }
  
  cat("\n=== 安装完成 ===\n")
  cat("成功:", success_count, "\n")
  cat("失败:", length(failed_list), "\n")
  if (length(failed_list) > 0) {
    cat("失败列表:\n", paste(failed_list, collapse = "\n"), "\n")
  }
}
 
# 使用示例
github_repos <- c("tpq/kpmt", "cBioPortal/cgdsr", "IOBR/IOBR")
batch_install_github_backup(github_repos)
 
cran_packages <- c("ggplot2", "dplyr", "limma")
batch_install_cran_backup(cran_packages)
