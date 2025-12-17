import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, BarycentricInterpolator

# ==========================================
# 配置中文显示
# ==========================================
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题

def main():
    # 1. 定义函数与节点
    # ------------------------------------------
    def f(x):
        return 1.0 / (1.0 + 25.0 * x**2)

    start, end = -5.0, 5.0
    n = 10  # 插值区间数 (节点数 = n+1)
    
    # 生成等距节点
    x_nodes = np.linspace(start, end, n + 1)
    y_nodes = f(x_nodes)
    
    # 生成用于绘图的密集点
    x_plot = np.linspace(start, end, 500)
    y_true = f(x_plot)

    # 2. 计算插值
    # ------------------------------------------
    # 拉格朗日插值 (使用重心插值法实现，数值更稳定)
    poly_lagrange = BarycentricInterpolator(x_nodes, y_nodes)
    y_lagrange = poly_lagrange(x_plot)
    
    # 三次样条插值
    cs = CubicSpline(x_nodes, y_nodes)
    y_spline = cs(x_plot)

    # 3. 绘图对比
    # ------------------------------------------
    plt.figure(figsize=(10, 6), dpi=120)
    
    # 真实函数
    plt.plot(x_plot, y_true, 'k--', linewidth=1.5, label='原函数 $f(x)$', alpha=0.6)
    
    # 拉格朗日插值
    plt.plot(x_plot, y_lagrange, 'r-', linewidth=1.5, label=f'拉格朗日插值 (n={n})')
    
    # 三次样条插值
    plt.plot(x_plot, y_spline, 'g-', linewidth=2, label=f'三次样条插值 (n={n})')
    
    # 插值节点
    plt.scatter(x_nodes, y_nodes, color='black', s=40, zorder=5, label='插值节点')

    # 装饰
    plt.title(f"龙格现象分析 (节点数 n={n})", fontsize=16)
    plt.xlabel('自变量 x', fontsize=12)
    plt.ylabel('函数值 y', fontsize=12)
    plt.ylim(-0.5, 1.2)  # 限制Y轴范围以突显振荡
    plt.legend(loc='upper right')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('runge_comparison_cn.png') # 保存为新文件名
    plt.show()

if __name__ == "__main__":
    main()