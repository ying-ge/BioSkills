#!/usr/bin/env python3
"""
Python包管理脚本 - 用于本地安装和管理下载的包
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import json

class PackageManager:
    def __init__(self, backup_dir="Python_packages_backup"):
        self.backup_dir = Path(backup_dir)
        
    def list_packages(self, category=None):
        """列出可用包"""
        if category:
            packages_dir = self.backup_dir / category
            if not packages_dir.exists():
                print(f"❌ Category not found: {category}")
                return
        else:
            packages_dir = self.backup_dir
            
        if not packages_dir.exists():
            print(f"❌ Backup directory not found: {self.backup_dir}")
            print("💡 Make sure you've run the download workflow first")
            return
            
        print("🚀 Available Python Packages")
        print("=" * 40)
        
        for category_dir in sorted(packages_dir.iterdir()):
            if category_dir.is_dir() and category_dir.name != "source":
                wheels = list(category_dir.glob("*.whl"))
                sources = list(category_dir.glob("source/*.tar.gz"))
                
                if wheels or sources:
                    print(f"\n📦 {category_dir.name}:")
                    
                    # 显示wheel包
                    for pkg in sorted(wheels):
                        info_file = category_dir / f"{pkg.stem.split('-')[0]}_info.json"
                        version = "unknown"
                        if info_file.exists():
                            try:
                                with open(info_file) as f:
                                    info = json.load(f)
                                    version = info.get('version', 'unknown')
                            except:
                                pass
                        print(f"  📄 {pkg.name} (v{version})")
                    
                    # 显示源码包
                    for src in sorted(sources):
                        print(f"  📁 source/{src.name}")
                        
    def install_package(self, package_name, category, use_source=False, force=False):
        """安装指定包"""
        category_dir = self.backup_dir / category
        
        if not category_dir.exists():
            print(f"❌ Category not found: {category}")
            return False
            
        if use_source:
            pkg_pattern = f"source/{package_name}.tar.gz"
            pkg_files = list(category_dir.glob(pkg_pattern))
            if not pkg_files:
                # 尝试其他命名方式
                pkg_files = list(category_dir.glob(f"source/{package_name}-*.tar.gz"))
        else:
            pkg_files = list(category_dir.glob(f"{package_name}-*.whl"))
            
        if not pkg_files:
            print(f"❌ Package not found: {package_name} in category {category}")
            return False
            
        pkg_path = pkg_files[0]  # 使用第一个匹配的包
        
        # 构建安装命令
        cmd = ["pip", "install"]
        if force:
            cmd.append("--force-reinstall")
        cmd.append(str(pkg_path))
        
        print(f"🚀 Installing: {package_name}")
        print(f"📄 From: {pkg_path}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode == 0:
                print(f"✅ Successfully installed: {package_name}")
                return True
            else:
                print(f"❌ Installation failed: {package_name}")
                print(f"Error: {result.stderr}")
                return False
        except Exception as e:
            print(f"❌ Error installing {package_name}: {e}")
            return False
            
    def install_category(self, category, use_source=False, force=False):
        """安装整个类别的包"""
        category_dir = self.backup_dir / category
        
        if not category_dir.exists():
            print(f"❌ Category not found: {category}")
            return
            
        if use_source:
            packages = list(category_dir.glob("source/*.tar.gz"))
        else:
            packages = list(category_dir.glob("*.whl"))
            
        if not packages:
            print(f"❌ No packages found in category: {category}")
            return
            
        print(f"🚀 Installing {len(packages)} packages from category: {category}")
        
        success_count = 0
        for pkg in sorted(packages):
            package_name = pkg.stem.split('-')[0]  # 提取包名
            if self.install_package(package_name, category, use_source, force):
                success_count += 1
                
        print(f"✅ Installed {success_count}/{len(packages)} packages")
        
    def get_package_info(self, package_name, category):
        """获取包信息"""
        category_dir = self.backup_dir / category
        info_file = category_dir / f"{package_name}_info.json"
        
        if not info_file.exists():
            print(f"❌ Package info not found: {package_name}")
            return
            
        try:
            with open(info_file) as f:
                info = json.load(f)
                
            print(f"📋 Package Information: {package_name}")
            print("=" * 30)
            print(f"Name: {info.get('name', 'N/A')}")
            print(f"Version: {info.get('version', 'N/A')}")
            print(f"Description: {info.get('description', 'N/A')}")
            print(f"Homepage: {info.get('homepage', 'N/A')}")
            print(f"Download Date: {info.get('download_date', 'N/A')}")
            
        except Exception as e:
            print(f"❌ Error reading package info: {e}")

def main():
    parser = argparse.ArgumentParser(description="Python包管理工具")
    parser.add_argument("--backup-dir", default="Python_packages_backup", 
                       help="Backup directory path (default: Python_packages_backup)")
    parser.add_argument("--list", action="store_true", help="列出可用包")
    parser.add_argument("--category", help="指定包类别")
    parser.add_argument("--install", help="安装指定包")
    parser.add_argument("--install-category", help="安装整个类别")
    parser.add_argument("--use-source", action="store_true", help="使用源码安装")
    parser.add_argument("--force", action="store_true", help="强制重新安装")
    parser.add_argument("--info", help="获取包信息")
    
    args = parser.parse_args()
    
    manager = PackageManager(args.backup_dir)
    
    if args.list:
        manager.list_packages(args.category)
    elif args.install:
        if not args.category:
            print("❌ Please specify --category when installing a specific package")
            sys.exit(1)
        manager.install_package(args.install, args.category, args.use_source, args.force)
    elif args.install_category:
        manager.install_category(args.install_category, args.use_source, args.force)
    elif args.info:
        if not args.category:
            print("❌ Please specify --category when getting package info")
            sys.exit(1)
        manager.get_package_info(args.info, args.category)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
