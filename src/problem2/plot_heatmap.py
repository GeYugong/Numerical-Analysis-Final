import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# 读取数据
csv_path = "data/temperature_field.csv"
if not os.path.exists(csv_path):
    print(f"Error: {csv_path} not found. Run the C++ solver first.")
    exit()

# 由于 Eigen 导出的 CSV 可能没有 header，我们直接读取
# 假设数据是 50x50 的矩阵
data = pd.read_csv(csv_path, header=None).values

# 设置网格坐标
Ny, Nx = data.shape
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

# 绘图
plt.figure(figsize=(10, 8))
# 填充等温图
cp = plt.contourf(X, Y, data, 50, cmap='inferno')
cbar = plt.colorbar(cp)
cbar.set_label('Temperature ($^\circ$C)', rotation=270, labelpad=15)

# 添加等温线 (Isolines)
levels = np.linspace(0, 100, 11)
cs = plt.contour(X, Y, data, levels=levels, colors='white', linewidths=0.8, alpha=0.7)
plt.clabel(cs, inline=True, fmt='%1.0f', fontsize=8)

plt.title(f'Steady State Temperature Distribution (FVM {Nx}x{Ny})')
plt.xlabel('x (m)')
plt.ylabel('y (m)')

# 保存图片
plt.savefig("docs/temperature_distribution.png", dpi=300)
print("Plot saved to docs/temperature_distribution.png")
plt.show()