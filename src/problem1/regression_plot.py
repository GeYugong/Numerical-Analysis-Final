import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 配置中文显示 (解决中文乱码和负号显示问题)
# ==========================================
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体
plt.rcParams['axes.unicode_minus'] = False    # 解决保存图像是负号'-'显示为方块的问题

def main():
    # 1. 数据准备
    # ------------------------------------------
    # 已知原始数据
    t_known = np.array([0, 1, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14], dtype=float)
    h_known = np.array([320.5, 374.8, 419.2, 455.6, 479.1, 494.5, 455.8, 419.5, 376.2, 319.8, 256.4, 178.9])

    # 2. 数据恢复 (拉格朗日插值)
    # ------------------------------------------
    def lagrange_interp(x_nodes, y_nodes, x_target):
        res = 0.0
        n = len(x_nodes)
        for i in range(n):
            term = y_nodes[i]
            for j in range(n):
                if i != j:
                    term *= (x_target - x_nodes[j]) / (x_nodes[i] - x_nodes[j])
            res += term
        return res

    # 恢复 t=6, 7 (利用 t=4, 5, 9)
    mask_g1 = np.isin(t_known, [4, 5, 9])
    h_6 = lagrange_interp(t_known[mask_g1], h_known[mask_g1], 6)
    h_7 = lagrange_interp(t_known[mask_g1], h_known[mask_g1], 7)

    # 恢复 t=8 (利用 t=5, 9, 10)
    mask_g2 = np.isin(t_known, [5, 9, 10])
    h_8 = lagrange_interp(t_known[mask_g2], h_known[mask_g2], 8)

    t_restored = [6, 7, 8]
    h_restored = [h_6, h_7, h_8]

    # 3. 回归模型与绘图
    # ------------------------------------------
    # 回归系数 (来自 C++ 程序输出)
    a, b, c = -5.00646, 60.1232, 319.697
    
    def model(t):
        return a * t**2 + b * t + c

    # 生成平滑曲线数据
    t_smooth = np.linspace(0, 17, 200)
    h_smooth = model(t_smooth)

    # 预测落地时间
    t_land = 16.000145

    # 绘图设置
    plt.figure(figsize=(10, 6), dpi=120)
    
    # 绘制拟合曲线
    plt.plot(t_smooth, h_smooth, 'r-', linewidth=2, label=f'拟合模型: $h(t)={a:.2f}t^2+{b:.2f}t+{c:.2f}$')
    
    # 绘制原始数据点
    plt.scatter(t_known, h_known, color='blue', s=50, label='原始数据 (Data)', zorder=3)
    
    # 绘制恢复的数据点
    plt.scatter(t_restored, h_restored, color='green', marker='x', s=80, linewidth=2, label='恢复数据 (插值补全)', zorder=3)
    
    # 标记落地时间
    plt.scatter([t_land], [0], color='red', marker='*', s=150, label=f'预测落地: t={t_land:.2f}s', zorder=4)

    # 装饰
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')
    
    plt.title('飞行器轨迹分析：数据恢复与二次回归', fontsize=16)
    plt.xlabel('时间 t (s)', fontsize=12)
    plt.ylabel('高度 h (m)', fontsize=12)
    plt.legend(loc='best') # 自动寻找最佳位置
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('regression_curve_cn.png') # 保存为新文件名
    plt.show()

if __name__ == "__main__":
    main()