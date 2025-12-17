# 数值分析课程大作业

**作者**: Ge Yugong  
**项目**: 数值分析综合实验 - 飞行器遥测数据分析与二维稳态热传导模拟

---

## 📚 项目概述

本项目包含两个独立的数值分析应用实例，涵盖插值、回归、方程求根、线性方程组求解、偏微分方程离散化等多个核心数值方法：

- **Problem 1**: 飞行器遥测数据分析与控制系统
- **Problem 2**: 基于有限体积法(FVM)的二维稳态热传导数值模拟

项目采用 C++ 实现核心数值算法，并使用 Python 进行数据可视化，展现了从数学建模、算法实现到工程应用的完整数值计算流程。

---

## 🗂️ 项目结构

```
Numerical-Analysis-Final/
│
├── data/                          # 数据文件
│   ├── flight_data.csv           # 飞行器遥测数据
│   └── flight_data.md            # 数据说明文档
│
├── docs/                          # 文档目录（可选）
│
├── include/                       # C++ 头文件库
│   ├── NumCpp.h                  # 基础数值算法库（无依赖）
│   └── NumUtils.h                # Eigen 封装库（Problem 2）
│
├── src/                          # 源代码目录
│   │
│   ├── problem1/                 # Problem 1: 飞行器遥测分析
│   │   ├── main.cpp             # 主程序
│   │   ├── regression_plot.py   # 回归分析可视化
│   │   ├── runge_phenomenon.py  # 龙格现象演示
│   │   ├── question(corrected).md
│   │   └── 飞行器遥测数据分析与控制系统.md
│   │
│   └── problem2/                 # Problem 2: 热传导模拟
│       ├── main.cpp             # FVM 求解器
│       ├── plot_result.py       # 温度场可视化
│       ├── temperature_field.csv # 求解结果（程序生成）
│       ├── question.md
│       └── 基于有限体积法(FVM)的二维稳态热传导数值模拟.md
│
├── .vscode/                      # VS Code 配置
│   ├── c_cpp_properties.json    # C++ IntelliSense 配置
│   └── tasks.json               # 编译任务配置
│
└── README.md                     # 本文件
```

---

## 🚀 Problem 1: 飞行器遥测数据分析与控制系统

### 问题背景

一架试验飞行器在返回阶段因强电磁干扰导致 $t=6, 7, 8$ 秒的高度数据丢失。需要完成以下任务：

1. **数据恢复**: 使用插值算法填补缺失数据
2. **运动建模**: 二次多项式回归建立运动方程 $h(t) = at^2 + bt + c$
3. **落地预测**: 牛顿法求解 $h(t) = 0$，预测落地时间
4. **推力解算**: SOR 迭代法求解线性方程组，计算引擎推力分配

### 核心算法

| 算法模块 | 实现方法 | 用途 |
|---------|---------|------|
| 数据插值 | Lagrange 插值 | 恢复 t=6,7,8 的高度数据 |
| 回归分析 | 最小二乘法（正规方程） | 拟合二次运动模型 |
| 非线性方程求根 | 牛顿迭代法 | 预测落地时间 |
| 线性方程组求解 | SOR 超松弛迭代 | 计算引擎推力 |
| 三次样条插值 | 自然边界条件 | 平滑数据处理 |

### 运行方式

```bash
# 编译
cd src/problem1
g++ -std=c++17 -I../../include main.cpp -o main.exe

# 运行
./main.exe

# 可视化
python regression_plot.py
python runge_phenomenon.py
```

### 关键结果

- **拟合方程**: $h(t) = -5.00646 t^2 + 60.1232 t + 319.697$
- **预测落地时间**: $t \approx 16.00$ 秒
- **引擎推力**: $x_1=3.0, x_2=4.0, x_3=-5.0$
- **MSE**: 0.624（拟合精度良好）

---

## 🔥 Problem 2: 二维稳态热传导数值模拟

### 问题背景

使用**有限体积法 (FVM)** 求解单位正方形区域 $[0,1] \times [0,1]$ 上的稳态热传导方程：

$$-\nabla^2 T = 0$$

**边界条件**:
- 左边界 (西): $T = 100°C$ (Dirichlet)
- 右边界 (东): $T = 0°C$ (Dirichlet)
- 上下边界 (北南): 绝热 $\frac{\partial T}{\partial n} = 0$ (Neumann)

### 数值方法

- **离散化方法**: 有限体积法 (Finite Volume Method)
- **网格**: $50 \times 50$ 均匀结构化网格
- **稀疏矩阵求解器**: Eigen 库的 `SimplicialLLT` (适用对称正定矩阵)
- **矩阵规模**: $2500 \times 2500$ 稀疏系统

### 技术特点

1. **守恒性保证**: FVM 基于积分守恒律，天然满足局部和全局能量守恒
2. **高效求解**: 使用 Eigen 稀疏矩阵存储（三元组初始化 + CSC 格式）
3. **内存优化**: 预分配 `triplets.reserve(N*5)` 避免动态扩容
4. **可视化**: Python matplotlib 生成高质量温度云图

### 运行方式

```bash
# 编译（需要 Eigen 库）
cd src/problem2
g++ -std=c++17 -I../../include -ID:/0software/eigen-5.0.0 main.cpp -o main.exe

# 运行
./main.exe

# 可视化
python plot_result.py
```

### 关键结果

- **求解时间**: 约 0.2 秒（2500 个未知数）
- **温度分布**: 热量从左边界（100°C）水平传导至右边界（0°C）
- **输出文件**: `temperature_field.csv` (50x50 矩阵)

---

## 🛠️ 技术栈

| 组件 | 技术 | 说明 |
|------|------|------|
| 编程语言 | C++17 | 核心数值算法实现 |
| 线性代数库 | Eigen 5.0.0 | 稀疏矩阵求解（Problem 2） |
| 可视化 | Python 3.8 + matplotlib | 绘制回归曲线、温度场 |
| 编译器 | MinGW g++ 8.1.0 | Windows 环境编译 |
| 开发环境 | VS Code | 配置文件见 `.vscode/` |

---

## 📦 依赖安装

### C++ 部分

1. **编译器**: MinGW-w64 或 MSVC
2. **Eigen 库** (仅 Problem 2 需要):
   - 下载: https://eigen.tuxfamily.org/
   - 解压到: `D:/0software/eigen-5.0.0` (或修改 `.vscode/tasks.json` 中的路径)
   - Eigen 是纯头文件库，无需编译

### Python 部分

```bash
pip install numpy matplotlib pandas
```

---

## 🔧 编译配置说明

### IntelliSense 配置

编辑 `.vscode/c_cpp_properties.json`，确保包含路径正确：

```json
{
  "includePath": [
    "${workspaceFolder}/**",
    "${workspaceFolder}/include",
    "D:/0software/eigen-5.0.0"  // Eigen 路径
  ],
  "compilerPath": "D:/0software/x86_64-8.1.0-release-posix-seh-rt_v6-rev0/mingw64/bin/g++.exe",
  "cppStandard": "c++17"
}
```

### 编译任务

使用 VS Code 内置任务（`Ctrl+Shift+B`）或手动编译：

```bash
# Problem 1 (无外部依赖)
g++ -std=c++17 -I./include src/problem1/main.cpp -o problem1.exe

# Problem 2 (需要 Eigen)
g++ -std=c++17 -I./include -I/path/to/eigen src/problem2/main.cpp -o problem2.exe
```

---

## 📊 核心算法实现

### NumCpp.h - 基础数值库

```cpp
class Interpolator {
    static double lagrange(...);      // 拉格朗日插值
    static double newton(...);         // 牛顿插值
    static double piecewiseLinear(...); // 分段线性插值
};

class CubicSpline {                    // 三次样条插值
    void solveMoments();
    double evaluate(double x);
};

class LinearAlgebra {
    static Vector solveGaussian(...);  // 高斯消元
    static Vector solveSOR(...);       // SOR 迭代
    static Matrix multiply(...);       // 矩阵运算
};

class Solvers {
    static double newtonMethod(...);   // 牛顿法求根
};
```

### NumUtils.h - Eigen 封装库

```cpp
namespace NumUtils {
    using SpMat = Eigen::SparseMatrix<double>;
    using Vector = Eigen::VectorXd;
    
    class LinearSolver {
        static bool solve(const SpMat& A, const Vector& b, Vector& x, 
                         SolverType type = SimplicialLLT);
    };
    
    void saveToCSV(const std::string& filename, const Matrix& matrix);
}
```

---

## 📈 数值方法对比

| 方法 | Problem 1 应用 | Problem 2 应用 |
|------|---------------|---------------|
| 插值法 | Lagrange (局部低次) | - |
| 回归分析 | 最小二乘法 | - |
| 线性求解 | SOR 迭代 (小规模) | Cholesky 分解 (大规模稀疏) |
| 非线性求解 | 牛顿法 | - |
| 偏微分方程 | - | 有限体积法 (FVM) |

---

## 🎓 学习要点

### Problem 1 重点

1. **插值选点策略**: 避免龙格现象，采用分段低次插值
2. **最小二乘法**: 通过正规方程 $A^T A \mathbf{x} = A^T \mathbf{b}$ 求解
3. **牛顿法 vs 二分法**: 导数信息加速收敛，适合外推
4. **SOR 参数调优**: 松弛因子 $\omega=1.25$ 加速收敛

### Problem 2 重点

1. **FVM 离散化**: 积分守恒律 → 代数方程组
2. **边界条件处理**: Dirichlet (修正源项) / Neumann (置零通量)
3. **稀疏矩阵存储**: 三元组 → CSC 格式，内存效率 99%+
4. **Eigen 求解器选择**: SPD 矩阵用 SimplicialLLT，速度最快

---

## 🐛 常见问题

### 1. "Eigen/Sparse: No such file or directory"

**解决方案**: 
- 确保 Eigen 已下载并解压
- 检查 `.vscode/c_cpp_properties.json` 和 `.vscode/tasks.json` 中的路径
- 重新加载 VS Code 窗口

### 2. Python 报错 "无法从源解析导入 pandas"

**解决方案**:
```bash
pip install pandas matplotlib numpy
```

### 3. 编译任务找不到头文件

**解决方案**:
- 使用绝对路径: `-ID:/0software/eigen-5.0.0`
- 避免使用 UNIX 风格路径: `-I/d/0software/...`

---

## 📝 开发日志

- **2025-12-17**: 项目初始化，完成 Problem 1 和 Problem 2 的算法实现
- **配置优化**: 解决 Eigen 路径问题，统一 Windows 编译环境
- **文档完善**: 添加详细的数学推导和实现报告

---

## 📖 参考资料

1. **数值分析**: 李庆扬《数值分析》第五版
2. **有限体积法**: H.K. Versteeg, "An Introduction to Computational Fluid Dynamics: The Finite Volume Method"
3. **Eigen 文档**: https://eigen.tuxfamily.org/dox/

---

## 👤 作者信息

- **姓名**: Ge Yugong
- **课程**: 数值分析
- **完成时间**: 2025 年 12 月

---

## 📜 许可证

本项目仅用于学习和教学目的。