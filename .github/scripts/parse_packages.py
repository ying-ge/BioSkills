#!/usr/bin/env python3
"""
解析包列表文件，按类别分组
"""

import yaml
import sys
from pathlib import Path

def parse_package_list(file_path):
    """解析包列表文件，按类别分组"""
    categories = {
        'singlecell': [],
        'machine_learning': [],
        'core': [],
        'visualization': [],
        'genomics': [],
        'immunology': [],
        'deep_learning': [],
        'data_processing': [],
        'uncategorized': []
    }
    
    # 关键词映射 - 优化：预编译为更快的查找方式
    category_keywords = {
        '单细胞': 'singlecell',
        '机器学习': 'machine_learning',
        '核心': 'core',
        '可视化': 'visualization',
        '基因组': 'genomics',
        '免疫': 'immunology',
        '深度学习': 'deep_learning',
        '数据处理': 'data_processing'
    }
    
    current_category = 'uncategorized'
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                # 跳过空行
                if not line:
                    continue
                
                # 检查是否是类别标题
                if line.startswith('#') and '包' in line:
                    # 优化：使用生成器表达式和next()代替循环
                    # 这样找到第一个匹配就立即停止，而不是检查所有关键词
                    current_category = next(
                        (category for keyword, category in category_keywords.items() if keyword in line),
                        'uncategorized'
                    )
                    print(f"Category found: {line} -> {current_category}")
                    continue
                
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                # 包名
                package_name = line.split()[0]
                if package_name:
                    categories[current_category].append(package_name)
                    print(f"Line {line_num}: {package_name} -> {current_category}")
                    
    except FileNotFoundError:
        print(f"Error: Package list file not found: {file_path}")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing package list: {e}")
        sys.exit(1)
    
    return categories

def main():
    if len(sys.argv) != 3:
        print("Usage: python parse_packages.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print(f"Parsing packages from: {input_file}")
    categories = parse_package_list(input_file)
    
    # 保存分类结果
    with open(output_file, 'w', encoding='utf-8') as f:
        yaml.dump(categories, f, default_flow_style=False, allow_unicode=True)
    
    print(f"✅ Package categorization completed")
    total_packages = sum(len(packages) for packages in categories.values())
    for category, packages in categories.items():
        print(f"📦 {category}: {len(packages)} packages")
    print(f"🎯 Total: {total_packages} packages")

if __name__ == "__main__":
    main()
