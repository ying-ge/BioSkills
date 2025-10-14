# scripts/extract_versions_from_figureya.R
#!/usr/bin/env Rscript
# 从FigureYa仓库提取包版本信息
 
suppressPackageStartupMessages({
  library(rvest)
  library(stringr)
  library(jsonlite)
  library(httr)
})
 
# 配置
FIGUREYA_REPO <- "ying-ge/FigureYa"
GITHUB_TOKEN <- Sys.getenv("GITHUB_TOKEN")
OUTPUT_FILE <- "metadata/figureya_packages.json"
 
# 从GitHub API获取文件列表
get_html_files_from_github <- function(repo, token = NULL) {
  headers <- list()
  if(!is.null(token) && token != "") {
    headers[["Authorization"]] <- paste("token", token)
  }
  
  # 递归获取所有HTML文件
  get_files_recursive <- function(path = "") {
    url <- paste0("https://api.github.com/repos/", repo, "/contents/", path)
    
    response <- GET(url, add_headers(.headers = headers))
    
    if(status_code(response) != 200) {
      cat("Failed to fetch:", url, "\n")
      return(character(0))
    }
    
    content <- content(response, "parsed")
    html_files <- character(0)
    
    for(item in content) {
      if(item$type == "file" && grepl("\\.html$", item$name)) {
        # 跳过docs目录
        if(!grepl("^docs/", item$path)) {
          html_files <- c(html_files, item$download_url)
        }
      } else if(item$type == "dir" && !item$name %in% c("docs", ".git", ".github")) {
        # 递归处理子目录
        sub_files <- get_files_recursive(item$path)
        html_files <- c(html_files, sub_files)
      }
    }
    
    return(html_files)
  }
  
  return(get_files_recursive())
}
 
# 从HTML URL提取session info
extract_session_info_from_url <- function(html_url) {
  tryCatch({
    response <- GET(html_url)
    if(status_code(response) != 200) return(NULL)
    
    html_content <- read_html(content(response, "text"))
    
    # 查找session info
    session_blocks <- html_content %>%
      html_nodes("pre") %>%
      html_text()
    
    session_text <- session_blocks[grepl("attached base packages|other attached packages", session_blocks)]
    
    if(length(session_text) == 0) return(NULL)
    
    packages <- list()
    
    # 解析R版本
    r_version <- NA
    r_version_match <- str_extract(session_text, "R version [0-9.]+")
    if(!is.na(r_version_match)) {
      r_version <- str_extract(r_version_match, "[0-9.]+")
    }
    
    # 提取attached packages
    if(grepl("other attached packages:", session_text)) {
      attached_section <- str_extract(session_text, "other attached packages:.*?(?=\\n\\n|loaded via|$)")
      attached_packages <- str_extract_all(attached_section, "\\w+_[0-9.-]+")[[1]]
      
      for(pkg_info in attached_packages) {
        parts <- str_split(pkg_info, "_")[[1]]
        if(length(parts) == 2) {
          packages[[parts[1]]] <- list(
            version = parts[2],
            source = "CRAN",
            priority = "attached",
            r_version = r_version
          )
        }
      }
    }
    
    # 提取loaded packages
    if(grepl("loaded via a namespace", session_text)) {
      loaded_section <- str_extract(session_text, "loaded via a namespace.*")
      loaded_packages <- str_extract_all(loaded_section, "\\w+_[0-9.-]+")[[1]]
      
      for(pkg_info in loaded_packages) {
        parts <- str_split(pkg_info, "_")[[1]]
        if(length(parts) == 2 && !parts[1] %in% names(packages)) {
          packages[[parts[1]]] <- list(
            version = parts[2],
            source = "CRAN",
            priority = "loaded",
            r_version = r_version
          )
        }
      }
    }
    
    return(list(
      packages = packages,
      r_version = r_version,
      source_url = html_url
    ))
    
  }, error = function(e) {
    cat("Error processing", html_url, ":", e$message, "\n")
    return(NULL)
  })
}
 
# 主函数
main <- function() {
  cat("从FigureYa仓库提取包版本信息...\n")
  
  # 获取所有HTML文件
  html_urls <- get_html_files_from_github(FIGUREYA_REPO, GITHUB_TOKEN)
  cat("找到", length(html_urls), "个HTML文件\n")
  
  # 收集所有包信息
  all_data <- list(
    packages = list(),
    metadata = list(
      extraction_date = Sys.time(),
      figureya_repo = FIGUREYA_REPO,
      total_html_files = length(html_urls),
      r_versions = c()
    )
  )
  
  for(i in seq_along(html_urls)) {
    cat("处理 (", i, "/", length(html_urls), "):", basename(html_urls[i]), "\n")
    
    session_data <- extract_session_info_from_url(html_urls[i])
    
    if(!is.null(session_data)) {
      # 记录R版本
      if(!is.na(session_data$r_version)) {
        all_data$metadata$r_versions <- unique(c(all_data$metadata$r_versions, session_data$r_version))
      }
      
      # 合并包信息
      for(pkg_name in names(session_data$packages)) {
        pkg_info <- session_data$packages[[pkg_name]]
        
        if(!pkg_name %in% names(all_data$packages)) {
          all_data$packages[[pkg_name]] <- list()
        }
        
        # 按R版本组织
        r_ver <- ifelse(is.na(pkg_info$r_version), "unknown", pkg_info$r_version)
        
        if(!r_ver %in% names(all_data$packages[[pkg_name]])) {
          all_data$packages[[pkg_name]][[r_ver]] <- pkg_info
          all_data$packages[[pkg_name]][[r_ver]]$source_files <- c()
        }
        
        all_data$packages[[pkg_name]][[r_ver]]$source_files <- 
          unique(c(all_data$packages[[pkg_name]][[r_ver]]$source_files, session_data$source_url))
      }
    }
    
    # 每10个文件保存一次进度
    if(i %% 10 == 0) {
      dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
      writeLines(toJSON(all_data, pretty = TRUE, auto_unbox = TRUE), 
                paste0(OUTPUT_FILE, ".tmp"))
    }
  }
  
  # 保存最终结果
  dir.create(dirname(OUTPUT_FILE), showWarnings = FALSE, recursive = TRUE)
  writeLines(toJSON(all_data, pretty = TRUE, auto_unbox = TRUE), OUTPUT_FILE)
  
  # 删除临时文件
  if(file.exists(paste0(OUTPUT_FILE, ".tmp"))) {
    file.remove(paste0(OUTPUT_FILE, ".tmp"))
  }
  
  cat("\n提取完成!\n")
  cat("唯一包数量:", length(all_data$packages), "\n")
  cat("R版本:", paste(all_data$metadata$r_versions, collapse = ", "), "\n")
}
 
if(!interactive()) {
  main()
}
