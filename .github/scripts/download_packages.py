#!/usr/bin/env python3
"""
下载Python包的wheel和源码
"""

import yaml
import subprocess
import os
import sys
import requests
import json
from pathlib import Path
import time

def download_package(package_name, category_dir):
    """下载指定包的wheel和源码"""
    try:
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
            
        # 获取包信息
        try:
            print(f"  ℹ️  Fetching info for {package_name}...")
            response = requests.get(f"https://pypi.org/pypi/{package_name}/json", timeout=10)
            if response.status_code == 200:
                pkg_info = response.json()
                version = pkg_info.get('info', {}).get('version', 'unknown')
                description = pkg_info.get('info', {}).get('summary', 'No description')
                homepage = pkg_info.get('info', {}).get('home_page', '')
                
                # 保存包信息
                info_file = wheel_dir / f"{package_name}_info.json"
                with open(info_file, 'w', encoding='utf-8') as f:
                    json.dump({
                        'name': package_name,
                        'version': version,
                        'description': description,
                        'homepage': homepage,
                        'download_date': os.environ.get('GITHUB_RUN_STARTED_AT', 'unknown')
                    }, f, indent=2, ensure_ascii=False)
                
                print(f"    ✅ Info saved: {package_name} v{version}")
            else:
                print(f"    ⚠️  Could not fetch info for {package_name}: HTTP {response.status_code}")
                
        except requests.exceptions.RequestException as e:
            print(f"    ⚠️  Network error fetching info for {package_name}: {e}")
        except Exception as e:
            print(f"    ⚠️  Error saving info for {package_name}: {e}")
            
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
    
    # 下载所有包
    total_packages = sum(len(packages) for packages in categories.values())
    downloaded_count = 0
    successful_downloads = 0
    
    for category, packages in categories.items():
        print(f"\n📦 Downloading {len(packages)} packages for category: {category}")
        category_dir = Path(output_dir) / category
        
        for i, package in enumerate(packages, 1):
            print(f"\n[{i}/{len(packages)}] 🔄 Processing: {package}")
            success_count = download_package(package, category_dir)
            downloaded_count += 1
            successful_downloads += success_count
            
            # 添加小延迟避免频繁请求
            time.sleep(0.5)
    
    print(f"\n🎉 Download completed!")
    print(f"📊 Statistics:")
    print(f"  Total packages processed: {downloaded_count}")
    print(f"  Successful downloads: {successful_downloads}")
    print(f"  Success rate: {successful_downloads/(downloaded_count*2)*100:.1f}%")

if __name__ == "__main__":
    main()
