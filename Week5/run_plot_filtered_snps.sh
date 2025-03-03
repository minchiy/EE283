#!/bin/bash
#SBATCH --job-name=plot_snp  # 設定作業名稱
#SBATCH --output=plot_snp.out  # 設定輸出檔案
#SBATCH --error=plot_snp.err  # 設定錯誤輸出檔案
#SBATCH --time=01:00:00  # 設定作業運行時間（這裡設定 1 小時）
#SBATCH --mem=4GB  # 設定記憶體使用限制

# 載入 Python 模組（視你的系統設定）
module load python/3.8

# 安裝必要的 Python 套件（如果需要）
pip install --user matplotlib numpy

# 執行 Python 腳本
python3 plot_filtered_snps.py
