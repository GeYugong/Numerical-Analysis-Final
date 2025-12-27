/**
 * File: src/problem1/main.cpp
 * Description: 飞行器数据恢复(Q1) + 龙格现象分析(附加题) 
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include "../../include/NumCpp.h"

using namespace std;

// 辅助函数：读取CSV文件 
vector<Point> readCSV(const string &filename)
{
    vector<Point> data;
    ifstream file(filename);
    string line;

    if (!file.is_open())
    {
        cerr << "错误：无法打开文件 " << filename << endl;
        return data;
    }

    // 跳过表头
    getline(file, line);

    while (getline(file, line))
    {
        stringstream ss(line);
        string t_str, h_str;

        if (getline(ss, t_str, ',') && getline(ss, h_str, ','))
        {
            double t = stod(t_str);
            try
            {
                double h = stod(h_str);
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
Point getPointByTime(const vector<Point> &all_data, double t)
{
    for (const auto &p : all_data)
    {
        if (abs(p.t - t) < 1e-5)
            return p;
    }
    return {t, 0.0};
}

// 多项式拟合函数
Vector performPolyFit(const vector<Point> &data, int degree)
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
            X[i][j] = pow(t, j);
        }
    }

    // 2. 构建正规方程
    Matrix XT = LinearAlgebra::transpose(X);
    Matrix A = LinearAlgebra::multiply(XT, X);
    Vector B = LinearAlgebra::multiply(XT, Y);

    // 3. 求解
    return LinearAlgebra::solveGaussian(A, B);
}


// 附加题逻辑: 龙格现象分析
void runBonusQuestion()
{
    cout << "\n" << endl;
    cout << "   正在运行附加题：龙格现象分析" << endl;
    cout << "\n" << endl;

    auto f = [](double x)
    { return 1.0 / (1.0 + 25.0 * x * x); };

    vector<int> n_values = {2, 4, 6, 8, 10};

    // 打印表头
    cout << left << setw(6) << "n" 
              << setw(15) << "Lagrange误差" 
              << setw(15) << "Newton误差" 
              << setw(15) << "分段线性误差" 
              << setw(15) << "三次样条误差" << endl;
    cout << string(70, '-') << endl;

    for (int n : n_values)
    {
        // 1. 生成插值节点
        vector<Point> nodes;
        double start = -5.0, end = 5.0;
        double step = (end - start) / n;

        for (int i = 0; i <= n; ++i)
        {
            double x = start + i * step;
            nodes.push_back({x, f(x)});
        }

        // 2. 准备三次样条对象
        CubicSpline spline(nodes);

        // 3. 计算误差 
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

            max_err_lagrange = max(max_err_lagrange, abs(y_lag - y_true));
            max_err_newton = max(max_err_newton, abs(y_new - y_true));
            max_err_piecewise = max(max_err_piecewise, abs(y_pw - y_true));
            max_err_spline = max(max_err_spline, abs(y_spl - y_true));
        }

        // 输出
        cout << left << setw(6) << n 
                  << setw(15) << max_err_lagrange 
                  << setw(15) << max_err_newton 
                  << setw(15) << max_err_piecewise 
                  << setw(15) << max_err_spline << endl;
    }
}

int main()
{
    // Part 1: 数据读取与插值恢复
    string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    vector<Point> flightData = readCSV(filepath);

    if (flightData.empty())
    {
        cerr << "加载数据错误。请检查 main.cpp 的文件路径。\n";
        return 1;
    }

    vector<Point> group1 = {getPointByTime(flightData, 4), getPointByTime(flightData, 5), getPointByTime(flightData, 9)};
    vector<Point> group2 = {getPointByTime(flightData, 5), getPointByTime(flightData, 9), getPointByTime(flightData, 10)};

    cout << "--- [步骤 1] 数据恢复 (插值算法) ---" << endl;
    
    double h6 = Interpolator::lagrange(group1, 6.0);
    double h7 = Interpolator::lagrange(group1, 7.0);
    double h8 = Interpolator::lagrange(group2, 8.0);

    flightData[6].h = h6;
    flightData[7].h = h7;
    flightData[8].h = h8;

    cout << "恢复 t=6: " << h6 << endl;
    cout << "恢复 t=7: " << h7 << endl;
    cout << "恢复 t=8: " << h8 << endl;
    cout << "数据恢复完成。\n";

    // Part 2: 模型构建与回归 (对比线性与二次模型)
    cout << "\n--- [步骤 2] 模型构建与回归 (多项式拟合对比) ---" << endl;
    
    // 1. 尝试线性拟合
    Vector coeffs_lin = performPolyFit(flightData, 1);
    double mse_lin = 0.0;
    for (const auto &p : flightData) {
        double pred = coeffs_lin[1] * p.t + coeffs_lin[0];
        mse_lin += pow(pred - p.h, 2);
    }
    mse_lin /= flightData.size();

    // 2. 尝试二次拟合
    Vector coeffs_quad = performPolyFit(flightData, 2);
    double mse_quad = 0.0;
    for (const auto &p : flightData) {
        double pred = coeffs_quad[2] * p.t * p.t + coeffs_quad[1] * p.t + coeffs_quad[0];
        mse_quad += pow(pred - p.h, 2);
    }
    mse_quad /= flightData.size();

    // 3. 输出对比结果与论证
    cout << "模型对比结果 (MSE):" << endl;
    cout << "  线性模型   MSE: " << mse_lin << endl;
    cout << "  二次模型   MSE: " << mse_quad << endl;

    cout << "\n论证结论：" << endl;
    if (mse_quad < mse_lin) {
        cout << "  数据显示二次模型的均方误差 (MSE) 显著小于线性模型。" << endl;
        cout << "  结合物理学原理 (h=at^2+bt+c)，二次模型更适合。" << endl;
    }

    double c = coeffs_quad[0];
    double b = coeffs_quad[1];
    double a = coeffs_quad[2];

    cout << "最终拟合函数表达式: h(t) = " << a << " * t^2 + " << b << " * t + " << c << endl;

    // Part 3: 落地时间预测 (Newton法)
    cout << "\n--- [步骤 3] 落地时间预测 ---" << endl;

    auto h_func = [&](double t) { return a * t * t + b * t + c; };
    auto h_deriv = [&](double t) { return 2 * a * t + b; };

    double landing_time = Solvers::newtonMethod(h_func, h_deriv, 14.0);

    cout << "预测落地时间 (h=0): " << setprecision(6) << fixed << landing_time << " 秒" << endl;

    // Part 4: 引擎推力解算 (SOR迭代法)
    cout << "\n--- [步骤 4] 引擎推力控制 (SOR 方法) ---" << endl;

    Matrix A_eng = {{4.0, 3.0, 0.0}, {3.0, 4.0, -1.0}, {0.0, -1.0, 4.0}};
    Vector b_eng = {24.0, 30.0, -24.0};

    pair<Vector, int> result = LinearAlgebra::solveSOR(A_eng, b_eng, 1.25, {1,1,1}, 1e-6, 100);
    Vector thrusts = result.first;
    int iterations = result.second;

    cout << "SOR 迭代收敛于 " << iterations << " 次迭代。" << endl;
    cout << "引擎推力解: x1=" << thrusts[0] << ", x2=" << thrusts[1] << ", x3=" << thrusts[2] << endl;

    // 运行附加题
    runBonusQuestion();
    return 0;
}