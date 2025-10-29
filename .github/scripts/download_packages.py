#!/usr/bin/env python3
"""
下载Python包的wheel和源码
优化版本：支持并发下载和缓存检查
"""

import yaml
import subprocess
import os
import sys
import requests
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache

@lru_cache(maxsize=1000)
def get_package_info(package_name):
    """获取包信息（带缓存）"""
    try:
        response = requests.get(f"https://pypi.org/pypi/{package_name}/json", timeout=10)
        if response.status_code == 200:
            pkg_info = response.json()
            return {
                'name': package_name,
                'version': pkg_info.get('info', {}).get('version', 'unknown'),
                'description': pkg_info.get('info', {}).get('summary', 'No description'),
                'homepage': pkg_info.get('info', {}).get('home_page', ''),
                'download_date': os.environ.get('GITHUB_RUN_STARTED_AT', 'unknown')
            }
    except Exception as e:
        print(f"    ⚠️  Error fetching info for {package_name}: {e}")
    return None


def check_package_exists(package_name, category_dir):
    """检查包是否已下载（避免重复下载）"""
    wheel_dir = Path(category_dir)
    source_dir = wheel_dir / "source"
    
    # 检查是否已有任何wheel或源码文件
    if wheel_dir.exists():
        # 检查wheel文件
        wheel_files = list(wheel_dir.glob(f"{package_name}*.whl"))
        source_files = []
        if source_dir.exists():
            source_files = list(source_dir.glob(f"{package_name}*"))
        
        if wheel_files or source_files:
            return True
    return False


def download_package(package_name, category_dir):
    """下载指定包的wheel和源码"""
    try:
        # 优化：检查是否已下载，避免重复下载
        if check_package_exists(package_name, category_dir):
            print(f"  ⏭️  Skipping {package_name} (already downloaded)")
            return 2  # 返回2表示已存在
        
        # 创建子目录
        wheel_dir = Path(category_dir)
        source_dir = Path(category_dir) / "source"
        wheel_dir.mkdir(parents=True, exist_ok=True)
        source_dir.mkdir(parents=True, exist_ok=True)
        
        success_count = 0
        
        # 下载wheel文件
        print(f"  📦 Downloading wheel for {package_name}...")
        cmd = f"pip download --no-deps {package_name} -d {wheel_dir}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"    ✅ Wheel downloaded: {package_name}")
            success_count += 1
        else:
            print(f"    ❌ Wheel download failed: {package_name}")
            print(f"    Error: {result.stderr}")
                
        # 下载源码
        print(f"  📄 Downloading source for {package_name}...")
        cmd = f"pip download --no-deps --no-binary :all: {package_name} -d {source_dir}"
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"    ✅ Source downloaded: {package_name}")
            success_count += 1
        else:
            print(f"    ❌ Source download failed: {package_name}")
            print(f"    Error: {result.stderr}")
            
        # 获取包信息（优化：使用缓存）
        pkg_info = get_package_info(package_name)
        if pkg_info:
            info_file = wheel_dir / f"{package_name}_info.json"
            with open(info_file, 'w', encoding='utf-8') as f:
                json.dump(pkg_info, f, indent=2, ensure_ascii=False)
            print(f"    ✅ Info saved: {package_name} v{pkg_info['version']}")
            
        return success_count
        
    except Exception as e:
        print(f"    ❌ Error downloading {package_name}: {e}")
        return 0

def main():
    if len(sys.argv) != 3:
        print("Usage: python download_packages.py <classified_packages.yml> <output_dir>")
        sys.exit(1)
    
    classified_file = sys.argv[1]
    output_dir = sys.argv[2]
    
    print(f"Downloading packages to: {output_dir}")
    
    # 读取分类后的包列表
    try:
        with open(classified_file, 'r', encoding='utf-8') as f:
            categories = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Error: Classified packages file not found: {classified_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading classified packages: {e}")
        sys.exit(1)
    
    # 优化：使用并发下载来加速I/O密集型操作
    # 使用ThreadPoolExecutor而不是顺序处理
    total_packages = sum(len(packages) for packages in categories.values())
    downloaded_count = 0
    successful_downloads = 0
    skipped_count = 0
    
    # 控制并发数量，避免过多并发导致的问题
    max_workers = min(10, os.cpu_count() * 2) if os.cpu_count() else 10
    
    print(f"📊 Total packages to process: {total_packages}")
    print(f"⚡ Using {max_workers} concurrent workers for downloads")
    
    for category, packages in categories.items():
        print(f"\n📦 Processing {len(packages)} packages for category: {category}")
        category_dir = Path(output_dir) / category
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有下载任务
            future_to_package = {
                executor.submit(download_package, package, category_dir): package
                for package in packages
            }
            
            # 处理完成的任务
            for i, future in enumerate(as_completed(future_to_package), 1):
                package = future_to_package[future]
                print(f"\n[{i}/{len(packages)}] 🔄 Processing: {package}")
                try:
                    success_count = future.result()
                    downloaded_count += 1
                    if success_count == 2:  # 已存在
                        skipped_count += 1
                    else:
                        successful_downloads += success_count
                except Exception as e:
                    print(f"    ❌ Exception processing {package}: {e}")
    
    print(f"\n🎉 Download completed!")
    print(f"📊 Statistics:")
    print(f"  Total packages processed: {downloaded_count}")
    print(f"  Packages skipped (already downloaded): {skipped_count}")
    print(f"  New successful downloads: {successful_downloads}")
    if downloaded_count > 0:
        print(f"  Success rate: {successful_downloads/((downloaded_count-skipped_count)*2)*100:.1f}%")

if __name__ == "__main__":
    main()
