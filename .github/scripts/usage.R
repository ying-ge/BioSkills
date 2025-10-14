# 完整的备份包使用示例
library(devtools)  # 如果没有需要先安装
 
# 1. 查看备份情况
list_github_backups()
list_cran_backups()
 
# 2. 安装单个包
install_from_github_backup("IOBR/IOBR")
install_from_cran_backup("ggplot2")
 
# 3. 批量安装
github_list <- c("tpq/kpmt", "zzwch/crosslink")
cran_list <- c("dplyr", "tidyr", "stringr")
 
batch_install_github_backup(github_list)
batch_install_cran_backup(cran_list)
 
# 4. 创建和使用本地仓库
create_local_repo()
use_local_repo()
 
# 现在可以像平常一样使用install.packages()
install.packages("ggplot2")  # 会优先从本地仓库安装
