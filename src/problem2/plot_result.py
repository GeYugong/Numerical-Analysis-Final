import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# ==========================================
# 核心修改：配置字体以支持中文显示
# ==========================================
# Windows 系统通常使用 'SimHei' (黑体) 或 'Microsoft YaHei' (微软雅黑)
plt.rcParams['font.sans-serif'] = ['SimHei'] 
# 解决负号 '-' 显示为方块的问题
plt.rcParams['axes.unicode_minus'] = False 

# 1. 读取数据 (使用自动路径，防止找不到文件)
current_dir = os.path.dirname(os.path.abspath(__file__))
csv_path = os.path.join(current_dir, "temperature_field.csv")

try:
    data = pd.read_csv(csv_path, header=None).values
except FileNotFoundError:
    print(f"错误：在 {csv_path} 找不到数据文件，请先运行 C++ 求解器。")
    exit()

# 2. 准备坐标网格
Ny, Nx = data.shape
x = np.linspace(0, 1, Nx)
y = np.linspace(0, 1, Ny)
X, Y = np.meshgrid(x, y)

# 3. 绘图
plt.figure(figsize=(10, 8))

# 绘制彩色云图
# cmap 依然推荐用 'inferno' 或 'jet'，视觉效果好
cp = plt.contourf(X, Y, data, 50, cmap='inferno') 
cbar = plt.colorbar(cp)
# 修改：Colorbar 标签改为中文
cbar.set_label('温度 ($^\circ$C)', rotation=270, labelpad=15, fontsize=12)

# 绘制等温线
cs = plt.contour(X, Y, data, 10, colors='white', alpha=0.5, linewidths=1)
plt.clabel(cs, inline=True, fontsize=10, fmt='%1.0f')

# 修改：标题和坐标轴标签改为中文
plt.title(f'稳态温度场分布 ({Nx}x{Ny} 均匀网格)', fontsize=14)
plt.xlabel('x (位置)', fontsize=12)
plt.ylabel('y (位置)', fontsize=12)

# 4. 保存并显示
save_path = os.path.join(current_dir, "result_plot_cn.png")
plt.savefig(save_path, dpi=300)
print(f"中文可视化结果已保存为: {save_path}")
plt.show()