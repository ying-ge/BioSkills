# 从文件读取GitHub仓库列表并备份
library(jsonlite)
 
# 读取仓库列表函数
read_github_repos <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("文件不存在: ", file_path)
  }
  
  # 读取文件内容
  lines <- readLines(file_path, warn = FALSE)
  
  # 过滤空行和注释行
  lines <- lines[lines != "" & !grepl("^#", lines)]
  
  # 提取GitHub仓库路径（支持完整URL或者简短格式）
  repos <- character(0)
  for (line in lines) {
    if (grepl("github.com", line)) {
      # 从完整URL提取仓库路径
      repo <- gsub(".*github.com/([^/]+/[^/]+).*", "\\1", line)
      repos <- c(repos, repo)
    } else if (grepl("^[^/]+/[^/]+$", line)) {
      # 直接是仓库路径格式
      repos <- c(repos, line)
    }
  }
  
  return(unique(repos))
}
 
# 增强版备份函数
backup_repos_from_file <- function(repo_file, backup_dir = "R_packages_backup") {
  # 读取仓库列表
  cat("正在从文件读取仓库列表:", repo_file, "\n")
  github_repos <- read_github_repos(repo_file)
  cat("找到", length(github_repos), "个仓库\n\n")
  
  # 创建备份目录
  if (!dir.exists(backup_dir)) {
    dir.create(backup_dir, recursive = TRUE)
  }
  
  backup_log <- data.frame(
    repo = character(),
    status = character(),
    download_date = character(),
    file_size = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (repo in github_repos) {
    repo_name <- gsub("/", "_", repo)
    repo_dir <- file.path(backup_dir, repo_name)
    if (!dir.exists(repo_dir)) {
      dir.create(repo_dir)
    }
    
    cat("正在备份:", repo, "\n")
    
    tryCatch({
      # 尝试下载主分支
      zip_urls <- c(
        paste0("https://github.com/", repo, "/archive/refs/heads/master.zip"),
        paste0("https://github.com/", repo, "/archive/refs/heads/main.zip")
      )
      
      downloaded <- FALSE
      final_file <- NULL
      
      for (i in seq_along(zip_urls)) {
        branch_name <- c("master", "main")[i]
        zip_file <- file.path(repo_dir, paste0(repo_name, "_", branch_name, ".zip"))
        
        tryCatch({
          download.file(zip_urls[i], zip_file, mode = "wb", quiet = TRUE)
          if (file.exists(zip_file) && file.size(zip_file) > 1000) {
            downloaded <- TRUE
            final_file <- zip_file
            break
          } else if (file.exists(zip_file)) {
            file.remove(zip_file)  # 删除失败的下载
          }
        }, error = function(e) NULL)
      }
      
      if (downloaded) {
        # 获取仓库基本信息
        api_url <- paste0("https://api.github.com/repos/", repo)
        repo_info <- tryCatch({
          fromJSON(api_url)
        }, error = function(e) NULL)
        
        if (!is.null(repo_info)) {
          meta_file <- file.path(repo_dir, "repository_info.json")
          writeLines(toJSON(repo_info, pretty = TRUE), meta_file)
        }
        
        backup_log <<- rbind(backup_log, data.frame(
          repo = repo,
          status = "success",
          download_date = as.character(Sys.Date()),
          file_size = file.size(final_file)
        ))
        
        cat("✓ 成功备份:", repo, "\n")
      } else {
        backup_log <<- rbind(backup_log, data.frame(
          repo = repo,
          status = "failed - no accessible branch",
          download_date = as.character(Sys.Date()),
          file_size = 0
        ))
        cat("✗ 备份失败:", repo, "\n")
      }
      
    }, error = function(e) {
      backup_log <<- rbind(backup_log, data.frame(
        repo = repo,
        status = paste("error:", e$message),
        download_date = as.character(Sys.Date()),
        file_size = 0
      ))
      cat("✗ 备份出错:", repo, "-", e$message, "\n")
    })
    
    # 避免GitHub API限制
    Sys.sleep(0.5)
  }
  
  # 保存备份日志
  log_file <- file.path(backup_dir, paste0("backup_log_", Sys.Date(), ".csv"))
  write.csv(backup_log, log_file, row.names = FALSE)
  
  # 生成备份报告
  success_count <- sum(backup_log$status == "success")
  cat("\n=== 备份完成 ===\n")
  cat("总数:", nrow(backup_log), "\n")
  cat("成功:", success_count, "\n")
  cat("失败:", nrow(backup_log) - success_count, "\n")
  cat("备份目录:", normalizePath(backup_dir), "\n")
  cat("日志文件:", log_file, "\n")
  
  return(backup_log)
}
 
# 使用方法
repo_file <- ".github/docs/github_library.txt"  # 或者你实际的文件路径
backup_log <- backup_repos_from_file(repo_file)
 
# 查看备份结果
print(backup_log)
