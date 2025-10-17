#!/usr/bin/env python3
"""
创建下载包的元数据和安装脚本
"""

import yaml
import json
import os
import sys
from pathlib import Path
from datetime import datetime

def create_metadata(classified_file, output_dir):
    """创建包元数据"""
    try:
        with open(classified_file, 'r') as f:
            categories = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading classified packages: {e}")
        return

    # 统计包信息
    metadata = {
        'download_info': {
            'date': datetime.utcnow().isoformat(),
            'python_version': os.environ.get('PYTHON_VERSION', 'unknown'),
            'trigger': os.environ.get('GITHUB_EVENT_NAME', 'manual'),
            'workflow_run_id': os.environ.get('GITHUB_RUN_ID', 'unknown'),
            'package_file': os.environ.get('PACKAGE_FILE', 'unknown'),
            'repository': os.environ.get('GITHUB_REPOSITORY', 'unknown')
        },
        'package_stats': {
            'total_packages': sum(len(packages) for packages in categories.values()),
            'categories': {cat: len(packages) for cat, packages in categories.items()}
        },
        'categories': {
            'singlecell': 'Single-cell RNA-seq analysis tools',
            'machine_learning': 'Machine learning and deep learning frameworks',
            'core': 'Core data science libraries',
            'visualization': 'Plotting and visualization tools',
            'genomics': 'Genomics and bioinformatics tools',
            'immunology': 'Immunology analysis packages',
            'deep_learning': 'Deep learning frameworks',
            'data_processing': 'Data processing and manipulation tools',
            'uncategorized': 'Other packages'
        }
    }

    # 保存metadata
    metadata_file = Path(output_dir) / 'metadata.yml'
    with open(metadata_file, 'w', encoding='utf-8') as f:
        yaml.dump(metadata, f, default_flow_style=False, allow_unicode=True)
    
    print(f"✅ Metadata saved to: {metadata_file}")
    return metadata

def create_install_script(output_dir, metadata):
    """创建安装助手脚本"""
    install_script = f'''#!/bin/bash
# Python包安装助手
# 生成时间: {datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')}
# 包总数: {metadata['package_stats']['total_packages']}
# Python版本: {metadata['download_info']['python_version']}

set -e

BACKUP_DIR="{output_dir}"

echo "🚀 Python包安装助手"
echo "===================="
echo "包总数: {metadata['package_stats']['total_packages']}"
echo "下载时间: {metadata['download_info']['date']}"
echo ""

# 显示可用类别
echo "📦 可用包类别:"
for category in "$BACKUP_DIR"/*; do
    if [ -d "$category" ]; then
        cat_name=$(basename "$category")
        if [[ "$cat_name" != "source" ]]; then
            package_count=$(find "$category" -name "*.whl" 2>/dev/null | wc -l)
            echo "  - $cat_name ($package_count packages)"
        fi
    fi
done

echo ""
echo "🎯 使用方法:"
echo "  # 安装特定类别的所有包"
echo "  pip install $BACKUP_DIR/singlecell/*.whl"
echo ""
echo "  # 安装特定包"
echo "  pip install $BACKUP_DIR/core/numpy*.whl"
echo ""
echo "  # 从源码安装"
echo "  pip install $BACKUP_DIR/core/source/numpy.tar.gz"
echo ""
echo "  # 查看包信息"
echo "  cat $BACKUP_DIR/core/numpy_info.json"
echo ""
echo "📋 按类别安装示例:"
'''

    # 添加类别安装示例
    for category, count in metadata['package_stats']['categories'].items():
        if count > 0:
            install_script += f'\n  # {metadata["categories"][category]} ({count} 包)\n'
            install_script += f'  # pip install $BACKUP_DIR/{category}/*.whl\n'

    install_script += '''
echo ""
echo "⚡ 快速安装命令:"
echo "  # 核心包"
echo "  pip install ''' + output_dir + '''/core/*.whl ''' + output_dir + '''/visualization/*.whl"
echo ""
echo "🔧 故障排除:"
echo "  # 如果wheel文件不兼容，使用源码安装"
echo "  pip install --no-binary :all: <package_name>"
echo ""
echo "📖 更多信息请查看: ''' + output_dir + '''/metadata.yml"
'''

    install_script_path = Path(output_dir) / 'install_packages.sh'
    with open(install_script_path, 'w') as f:
        f.write(install_script)
    
    # 设置可执行权限
    os.chmod(install_script_path, 0o755)
    print(f"✅ Install script created: {install_script_path}")

def create_report(classified_file, output_dir, metadata):
    """创建下载报告"""
    report_path = Path(output_dir) / 'download_report.md'
    
    try:
        with open(classified_file, 'r') as f:
            categories = yaml.safe_load(f)
    except Exception as e:
        print(f"Error reading classified packages for report: {e}")
        return

    report_content = f'''# Python Packages Download Report

## 📥 Download Information
- **Date**: {metadata['download_info']['date']}
- **Python Version**: {metadata['download_info']['python_version']}
- **Trigger**: {metadata['download_info']['trigger']}
- **Package File**: {metadata['download_info']['package_file']}
- **Workflow Run**: {metadata['download_info']['workflow_run_id']}
- **Repository**: {metadata['download_info']['repository']}

## 📊 Package Statistics
- **Total Packages**: {metadata['package_stats']['total_packages']}

### By Category:
| Category | Package Count | Examples |
|----------|---------------|----------|
'''

    for category, packages in categories.items():
        examples = ', '.join(packages[:3]) + ('...' if len(packages) > 3 else '')
        report_content += f'| {category} | {len(packages)} | {examples} |\n'

    report_content += f'''
## 🚀 Quick Start

### Install from downloaded packages:
```bash
# Install all single-cell packages
pip install {output_dir}/singlecell/*.whl

# Install specific package
pip install {output_dir}/core/scanpy*.whl

# Install from source
pip install {output_dir}/singlecell/source/scanpy.tar.gz
```

### Use installation helper:
```bash
./{output_dir}/install_packages.sh
```

## 📁 Directory Structure
```
{output_dir}/
├── singlecell/           # Single-cell analysis tools
│   ├── *.whl            # Wheel files
│   ├── source/          # Source archives
│   └── *_info.json      # Package information
├── machine_learning/     # ML frameworks
├── core/                # Core libraries
├── visualization/        # Plotting tools
├── genomics/           # Bioinformatics tools
├── immunology/         # Immunology packages
├── deep_learning/      # Deep learning frameworks
├── data_processing/    # Data processing tools
├── metadata.yml        # Package metadata
├── install_packages.sh # Installation helper
└── download_report.md  # This report
```

## 🛠️ Package Management Tools

### Local Management Script:
```python
# Copy the management script to your local machine
python .github/scripts/manage_packages.py --list
python .github/scripts/manage_packages.py --install scanpy --category singlecell
```

### Version Information:
All packages include version information in `*_info.json` files for reproducibility.

---
*Generated by GitHub Actions on {datetime.utcnow().strftime('%Y-%m-%d')}*
'''

    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report_content)
    
    print(f"✅ Report created: {report_path}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python create_metadata.py <classified_packages.yml> <output_dir>")
        sys.exit(1)
    
    classified_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    print("Creating metadata and reports...")
    
    metadata = create_metadata(classified_file, output_dir)
    create_install_script(output_dir, metadata)
    create_report(classified_file, output_dir, metadata)
    
    print("✅ All metadata and scripts created successfully!")

if __name__ == "__main__":
    main()
