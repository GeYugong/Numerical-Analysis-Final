/**
 * File: src/problem1/main.cpp
 * Description: 飞行器数据恢复(Q1) + 龙格现象分析(附加题) - 纯控制台输出版
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
// Windows环境下为了正确显示中文可能需要包含此头文件
// #include <windows.h> 
#include "../../include/NumCpp.h"

// 辅助函数：读取CSV文件 (这个必须保留，因为要读取输入数据)
std::vector<Point> readCSV(const std::string &filename)
{
    std::vector<Point> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "错误：无法打开文件 " << filename << std::endl;
        return data;
    }

    // 跳过表头
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string t_str, h_str;

        if (std::getline(ss, t_str, ',') && std::getline(ss, h_str, ','))
        {
            double t = std::stod(t_str);
            try
            {
                double h = std::stod(h_str);
                data.push_back({t, h});
            }
            catch (...)
            {
                data.push_back({t, -9999.0});
            }
        }
    }
    return data;
}

// 辅助函数：通过时间t查找高度h
Point getPointByTime(const std::vector<Point> &all_data, double t)
{
    for (const auto &p : all_data)
    {
        if (std::abs(p.t - t) < 1e-5)
            return p;
    }
    return {t, 0.0};
}

// 多项式拟合函数
Vector performPolyFit(const std::vector<Point> &data, int degree)
{
    int n = data.size();
    int terms = degree + 1;

    // 1. 构建设计矩阵 X
    Matrix X = LinearAlgebra::zeros(n, terms);
    Vector Y(n);

    for (int i = 0; i < n; ++i)
    {
        double t = data[i].t;
        Y[i] = data[i].h;
        for (int j = 0; j < terms; ++j)
        {
            X[i][j] = std::pow(t, j);
        }
    }

    // 2. 构建正规方程
    Matrix XT = LinearAlgebra::transpose(X);
    Matrix A = LinearAlgebra::multiply(XT, X);
    Vector B = LinearAlgebra::multiply(XT, Y);

    // 3. 求解
    return LinearAlgebra::solveGaussian(A, B);
}

// ==========================================
// 附加题逻辑: 龙格现象分析
// ==========================================
void runBonusQuestion()
{
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   正在运行附加题：龙格现象分析 (Runge's Phenomenon)" << std::endl;
    std::cout << "==========================================" << std::endl;

    auto f = [](double x)
    { return 1.0 / (1.0 + 25.0 * x * x); };

    std::vector<int> n_values = {2, 4, 6, 8, 10};

    // 打印表头
    std::cout << std::left << std::setw(6) << "n" 
              << std::setw(15) << "Lagrange误差" 
              << std::setw(15) << "Newton误差" 
              << std::setw(15) << "分段线性误差" 
              << std::setw(15) << "三次样条误差" << std::endl;
    std::cout << std::string(70, '-') << std::endl;

    for (int n : n_values)
    {
        // 1. 生成插值节点
        std::vector<Point> nodes;
        double start = -5.0, end = 5.0;
        double step = (end - start) / n;

        for (int i = 0; i <= n; ++i)
        {
            double x = start + i * step;
            nodes.push_back({x, f(x)});
        }

        // 2. 准备三次样条对象
        CubicSpline spline(nodes);

        // 3. 计算误差 (不再写入文件，仅计算最大误差)
        double plot_step = (end - start) / 200.0;
        double max_err_lagrange = 0.0;
        double max_err_newton = 0.0;
        double max_err_piecewise = 0.0;
        double max_err_spline = 0.0;

        for (int i = 0; i <= 200; ++i)
        {
            double x = start + i * plot_step;
            double y_true = f(x);

            double y_lag = Interpolator::lagrange(nodes, x);
            double y_new = Interpolator::newton(nodes, x);
            double y_pw = Interpolator::piecewiseLinear(nodes, x);
            double y_spl = spline.evaluate(x);

            max_err_lagrange = std::max(max_err_lagrange, std::abs(y_lag - y_true));
            max_err_newton = std::max(max_err_newton, std::abs(y_new - y_true));
            max_err_piecewise = std::max(max_err_piecewise, std::abs(y_pw - y_true));
            max_err_spline = std::max(max_err_spline, std::abs(y_spl - y_true));
        }

        // 表格形式输出
        std::cout << std::left << std::setw(6) << n 
                  << std::setw(15) << max_err_lagrange 
                  << std::setw(15) << max_err_newton 
                  << std::setw(15) << max_err_piecewise 
                  << std::setw(15) << max_err_spline << std::endl;
    }
}

int main()
{
    // ==========================================
    // Part 1: 数据读取与插值恢复
    // ==========================================
    // TODO: 请确保此处的路径与你的实际路径一致
    std::string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    std::vector<Point> flightData = readCSV(filepath);

    if (flightData.empty())
    {
        std::cerr << "加载数据错误。请检查 main.cpp 的文件路径。\n";
        return 1;
    }

    std::vector<Point> group1 = {getPointByTime(flightData, 4), getPointByTime(flightData, 5), getPointByTime(flightData, 9)};
    std::vector<Point> group2 = {getPointByTime(flightData, 5), getPointByTime(flightData, 9), getPointByTime(flightData, 10)};

    std::cout << "--- [步骤 1] 数据恢复 (插值算法) ---" << std::endl;
    
    double h6 = Interpolator::lagrange(group1, 6.0);
    double h7 = Interpolator::lagrange(group1, 7.0);
    double h8 = Interpolator::lagrange(group2, 8.0);

    flightData[6].h = h6;
    flightData[7].h = h7;
    flightData[8].h = h8;

    std::cout << "恢复 t=6: " << h6 << std::endl;
    std::cout << "恢复 t=7: " << h7 << std::endl;
    std::cout << "恢复 t=8: " << h8 << std::endl;
    std::cout << "数据恢复完成。\n";

    // ==========================================
    // Part 2: 模型构建与回归 (对比线性与二次模型)
    // ==========================================
    std::cout << "\n--- [步骤 2] 模型构建与回归 (多项式拟合对比) ---" << std::endl;
    
    // 1. 尝试线性拟合
    Vector coeffs_lin = performPolyFit(flightData, 1);
    double mse_lin = 0.0;
    for (const auto &p : flightData) {
        double pred = coeffs_lin[1] * p.t + coeffs_lin[0];
        mse_lin += std::pow(pred - p.h, 2);
    }
    mse_lin /= flightData.size();

    // 2. 尝试二次拟合
    Vector coeffs_quad = performPolyFit(flightData, 2);
    double mse_quad = 0.0;
    for (const auto &p : flightData) {
        double pred = coeffs_quad[2] * p.t * p.t + coeffs_quad[1] * p.t + coeffs_quad[0];
        mse_quad += std::pow(pred - p.h, 2);
    }
    mse_quad /= flightData.size();

    // 3. 输出对比结果与论证
    std::cout << "模型对比结果 (MSE):" << std::endl;
    std::cout << "  线性模型   MSE: " << mse_lin << std::endl;
    std::cout << "  二次模型   MSE: " << mse_quad << std::endl;

    std::cout << "\n论证结论：" << std::endl;
    if (mse_quad < mse_lin) {
        std::cout << "  数据显示二次模型的均方误差 (MSE) 显著小于线性模型。" << std::endl;
        std::cout << "  结合物理学原理 (h=at^2+bt+c)，二次模型更适合。" << std::endl;
    }

    double c = coeffs_quad[0];
    double b = coeffs_quad[1];
    double a = coeffs_quad[2];

    std::cout << "最终拟合函数表达式: h(t) = " << a << " * t^2 + " << b << " * t + " << c << std::endl;

    // ==========================================
    // Part 3: 落地时间预测 (Newton法)
    // ==========================================
    std::cout << "\n--- [步骤 3] 落地时间预测 ---" << std::endl;

    auto h_func = [&](double t) { return a * t * t + b * t + c; };
    auto h_deriv = [&](double t) { return 2 * a * t + b; };

    double landing_time = Solvers::newtonMethod(h_func, h_deriv, 14.0);

    std::cout << "预测落地时间 (h=0): " << std::setprecision(6) << std::fixed << landing_time << " 秒" << std::endl;

    // ==========================================
    // Part 4: 引擎推力解算 (SOR迭代法)
    // ==========================================
    std::cout << "\n--- [步骤 4] 引擎推力控制 (SOR 方法) ---" << std::endl;

    Matrix A_eng = {{4.0, 3.0, 0.0}, {3.0, 4.0, -1.0}, {0.0, -1.0, 4.0}};
    Vector b_eng = {24.0, 30.0, -24.0};

    std::pair<Vector, int> result = LinearAlgebra::solveSOR(A_eng, b_eng, 1.25, {1,1,1}, 1e-6, 100);
    Vector thrusts = result.first;
    int iterations = result.second;

    std::cout << "SOR 迭代收敛于 " << iterations << " 次迭代。" << std::endl;
    std::cout << "引擎推力解: x1=" << thrusts[0] << ", x2=" << thrusts[1] << ", x3=" << thrusts[2] << std::endl;

    // ------------------------------------------
    // 运行附加题
    // ------------------------------------------
    runBonusQuestion();
    return 0;
}