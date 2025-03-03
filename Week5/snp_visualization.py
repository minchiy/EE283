import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.colors import ListedColormap

# 設定檔案路徑
vcf_dir = "/pub/minchiy/EE283/Week5/SNPs"
genotype_file = f"{vcf_dir}/X_1Mb.012"
pos_file = f"{vcf_dir}/X_1Mb.012.pos"

# 檢查檔案是否存在
if not os.path.exists(genotype_file) or not os.path.exists(pos_file):
    print("Error: One or more input files are missing.")
    exit(1)

# 讀取基因型數據
geno_data = np.loadtxt(genotype_file, dtype=int)[:, 1:]  # 移除第一列索引

# 過濾缺失值（確保所有數據為 0, 1, 2）
geno_data[geno_data < 0] = 0

# 讀取 SNP 位點
snp_positions = np.loadtxt(pos_file, dtype=int)

# 設定顏色對應（0=綠, 1=黃, 2=紅）
cmap = ListedColormap(["green", "yellow", "red"])

# 繪製 SNP 可視化圖
plt.figure(figsize=(15, 5))
plt.imshow(geno_data, cmap=cmap, aspect="auto", interpolation="nearest")

# 添加 SNP 位置標籤
num_ticks = 10
tick_indices = np.linspace(0, len(snp_positions) - 1, num_ticks, dtype=int)
tick_labels = snp_positions[tick_indices]

plt.xticks(tick_indices, tick_labels, rotation=45)
plt.colorbar(label="Genotype (0=Ref/Ref, 1=Ref/Alt, 2=Alt/Alt)")
plt.xlabel("SNP Position on X Chromosome (1Mb)")
plt.ylabel("Samples")
plt.title("SNP Genotypes in X Chromosome (First 1Mb)")

# 存檔
output_path = f"{vcf_dir}/snp_visualization.png"
plt.savefig(output_path, dpi=300)
plt.show()

print(f"SNP visualization saved to {output_path}")

