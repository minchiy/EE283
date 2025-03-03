
import numpy as np
import matplotlib.pyplot as plt

# 設定檔案路徑
vcf_dir = "/pub/minchiy/EE283/Week5/SNPs"
genotype_file = f"{vcf_dir}/X_1Mb_filtered.012"

# 讀取基因型數據
geno_data = np.loadtxt(genotype_file, dtype=int)[:, 1:]  # 移除第一列索引

# 設定顏色對應（0=綠, 1=黃, 2=紅）
color_map = {0: "green", 1: "yellow", 2: "red"}
color_matrix = np.vectorize(color_map.get)(geno_data)

# 繪製 SNP 可視化圖
plt.figure(figsize=(15, 5))
plt.imshow(geno_data, cmap="viridis", aspect="auto", interpolation="nearest")
plt.colorbar(label="Genotype (0=Ref/Ref, 1=Ref/Alt, 2=Alt/Alt)")
plt.xlabel("SNP Index")
plt.ylabel("Samples")
plt.title("Filtered SNP Genotypes in X Chromosome (First 1Mb)")
plt.savefig(f"{vcf_dir}/filtered_snp_visualization.png", dpi=300)
plt.show()

