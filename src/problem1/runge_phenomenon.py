import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, BarycentricInterpolator, KroghInterpolator

# ==========================================
# 全局配置
# ==========================================
plt.rcParams['font.sans-serif'] = ['SimHei']  # 指定默认字体为黑体，解决中文乱码
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示为方块的问题

def main():
    # 1. 定义目标函数
    # ------------------------------------------
    def f(x):
        return 1.0 / (1.0 + 25.0 * x**2)

    start, end = -5.0, 5.0
    n_values = [2, 4, 6, 8, 10]  # 测试不同的节点数 (等分点数 = n+1)

    # 创建子图，共享X轴
    fig, axes = plt.subplots(len(n_values), 1, figsize=(10, 14), dpi=120, sharex=True)
    
    # 绘图用的密集点
    x_plot = np.linspace(start, end, 500)
    y_true = f(x_plot)

    legend_handles = None  # 用于存储图例句柄，保证只显示一次图例

    print(f"{'n':<4} | {'Lagrange':<12} | {'Newton':<12} | {'Piecewise':<12} | {'Spline':<12}")
    print("-" * 65)

    for ax, n in zip(axes, n_values):
        # 生成插值节点
        x_nodes = np.linspace(start, end, n + 1)
        y_nodes = f(x_nodes)

        # 构建插值模型
        # 注: 理论上 Lagrange 和 Newton 生成的是同一个多项式，曲线应重合
        lag_poly = BarycentricInterpolator(x_nodes, y_nodes)    # 重心插值 (对应 Lagrange)
        newton_poly = KroghInterpolator(x_nodes, y_nodes)       # Krogh插值 (对应 Newton形式)
        y_piecewise = np.interp(x_plot, x_nodes, y_nodes)       # 分段线性
        spline_obj = CubicSpline(x_nodes, y_nodes)              # 三次样条
        y_spline = spline_obj(x_plot)

        # 计算插值结果
        y_lagrange = lag_poly(x_plot)
        y_newton = newton_poly(x_plot)

        # 计算最大误差
        max_err_l = np.max(np.abs(y_lagrange - y_true))
        max_err_n = np.max(np.abs(y_newton - y_true))
        max_err_p = np.max(np.abs(y_piecewise - y_true))
        max_err_s = np.max(np.abs(y_spline - y_true))

        # 绘图
        line_true, = ax.plot(x_plot, y_true, 'k--', linewidth=1.5, alpha=0.6, label='原函数 (True)')
        line_lag,  = ax.plot(x_plot, y_lagrange, 'r-', linewidth=1.5, alpha=0.8, label='Lagrange')
        # Newton 线用点划线，防止完全覆盖 Lagrange (因为它们数学上是一样的)
        line_new,  = ax.plot(x_plot, y_newton, 'm-.', linewidth=1.5, alpha=0.8, label='Newton') 
        line_pw,   = ax.plot(x_plot, y_piecewise, 'b:', linewidth=1.5, label='分段线性 (Piecewise)')
        line_sp,   = ax.plot(x_plot, y_spline, 'g-', linewidth=2.0, label='三次样条 (Spline)')
        
        ax.scatter(x_nodes, y_nodes, color='black', s=40, zorder=5, label='插值节点')

        # 设置坐标轴范围 (限制Y轴以展示震荡)
        ax.set_ylim(-0.5, 1.2)
        ax.set_ylabel(f'n={n}', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        # 子标题显示误差
        ax.set_title(
            f'n={n} | MaxErr -> Lag: {max_err_l:.3f}, New: {max_err_n:.3f}, Pie: {max_err_p:.3f}, Spl: {max_err_s:.3f}',
            fontsize=10
        )

        if legend_handles is None:
            legend_handles = [line_true, line_lag, line_new, line_pw, line_sp]

        # 打印控制台信息
        print(f"{n:<4} | {max_err_l:<12.6f} | {max_err_n:<12.6f} | {max_err_p:<12.6f} | {max_err_s:<12.6f}")

    # 整体布局设置
    axes[-1].set_xlabel('自变量 x', fontsize=12)
    fig.suptitle("龙格现象分析：四种插值算法对比", fontsize=16, y=0.99)
    
    # 统一图例
    fig.legend(handles=legend_handles, loc='upper center', ncol=5, bbox_to_anchor=(0.5, 0.96))
    
    plt.tight_layout(rect=[0, 0, 1, 0.94]) # 留出顶部给标题和图例
    plt.savefig('runge_comparison_cn.png')
    plt.show()

if __name__ == "__main__":
    main()