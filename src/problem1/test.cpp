/**
 * File: src/problem1/main.cpp
 * Description: 飞行器数据恢复(Q1) + 龙格现象分析(Bonus)
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

// --- 辅助函数 (保持不变) ---
std::vector<Point> readCSV(const std::string &filename)
{
    std::vector<Point> data;
    std::ifstream file(filename);
    std::string line;
    if (!file.is_open())
        return data;
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
                data.push_back({t, std::stod(h_str)});
            }
            catch (...)
            {
                data.push_back({t, -9999.0});
            }
        }
    }
    return data;
}

Point getPointByTime(const std::vector<Point> &all_data, double t)
{
    for (const auto &p : all_data)
        if (std::abs(p.t - t) < 1e-5)
            return p;
    return {t, 0.0};
}

Vector performPolyFit(const std::vector<Point> &data, int degree)
{
    int n = data.size();
    int terms = degree + 1;
    Matrix X = LinearAlgebra::zeros(n, terms);
    Vector Y(n);
    for (int i = 0; i < n; ++i)
    {
        Y[i] = data[i].h;
        for (int j = 0; j < terms; ++j)
            X[i][j] = std::pow(data[i].t, j);
    }
    Matrix XT = LinearAlgebra::transpose(X);
    Matrix A = LinearAlgebra::multiply(XT, X);
    Vector B = LinearAlgebra::multiply(XT, Y);
    return LinearAlgebra::solveGaussian(A, B);
}

// ==========================================
// 附加题逻辑: 龙格现象分析
// ==========================================
void runBonusQuestion()
{
    std::cout << "\n==========================================" << std::endl;
    std::cout << "   Running Bonus Question: Runge's Phenomenon" << std::endl;
    std::cout << "==========================================" << std::endl;

    // 目标函数
    auto f = [](double x)
    { return 1.0 / (1.0 + 25.0 * x * x); };

    // 需要测试的节点数量 n (等分点数量 = n+1)
    std::vector<int> n_values = {2, 4, 6, 8, 10};

    for (int n : n_values)
    {
        // 1. 生成插值节点 (在 [-5, 5] 上等距分布)
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

        // 3. 在密集网格上计算插值结果并导出
        std::string filename = "data/bonus_runge_n" + std::to_string(n) + ".csv";
        std::ofstream outfile(filename);
        outfile << "x,True,Lagrange,Newton,Piecewise,Spline\n";

        // 生成 200 个绘图点
        double plot_step = (end - start) / 200.0;
        double max_err_lagrange = 0.0;
        double max_err_spline = 0.0;

        for (int i = 0; i <= 200; ++i)
        {
            double x = start + i * plot_step;
            double y_true = f(x);

            double y_lag = Interpolator::lagrange(nodes, x);
            double y_new = Interpolator::newton(nodes, x);
            double y_pw = Interpolator::piecewiseLinear(nodes, x);
            double y_spl = spline.evaluate(x);

            // 记录误差用于简单的控制台输出
            max_err_lagrange = std::max(max_err_lagrange, std::abs(y_lag - y_true));
            max_err_spline = std::max(max_err_spline, std::abs(y_spl - y_true));

            outfile << x << "," << y_true << "," << y_lag << "," << y_new << "," << y_pw << "," << y_spl << "\n";
        }

        std::cout << "n=" << n << ": Saved to " << filename << std::endl;
        std::cout << "     Max Error (Lagrange): " << max_err_lagrange << std::endl;
        std::cout << "     Max Error (Spline)  : " << max_err_spline << std::endl;
    }
}

int main()
{
    // ------------------------------------------
    // 第一题主任务 (保持原有代码，确保不破坏功能)
    // ------------------------------------------
    std::cout << "--- Main Task: Flight Data Analysis ---" << std::endl;
    std::string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    std::vector<Point> flightData = readCSV(filepath);

    if (!flightData.empty())
    {
        std::vector<Point> group1 = {getPointByTime(flightData, 4), getPointByTime(flightData, 5), getPointByTime(flightData, 9)};
        std::vector<Point> group2 = {getPointByTime(flightData, 5), getPointByTime(flightData, 9), getPointByTime(flightData, 10)};
        flightData[6].h = Interpolator::lagrange(group1, 6.0);
        flightData[7].h = Interpolator::lagrange(group1, 7.0);
        flightData[8].h = Interpolator::lagrange(group2, 8.0);

        Vector coeffs = performPolyFit(flightData, 2);
        auto h_func = [&](double t)
        { return coeffs[2] * t * t + coeffs[1] * t + coeffs[0]; };
        auto h_deriv = [&](double t)
        { return 2 * coeffs[2] * t + coeffs[1]; };
        double landing_time = Solvers::newtonMethod(h_func, h_deriv, 14.0);
        std::cout << "Landing Time: " << landing_time << " s" << std::endl;

        Matrix A_eng = {{4, 3, 0}, {3, 4, -1}, {0, -1, 4}};
        Vector b_eng = {24, 30, -24};
        auto res = LinearAlgebra::solveSOR(A_eng, b_eng, 1.25, {1, 1, 1}, 1e-6, 100);
        std::cout << "Engine: " << res.first[0] << ", " << res.first[1] << ", " << res.first[2] << std::endl;
    }

    // ------------------------------------------
    // 运行附加题
    // ------------------------------------------
    runBonusQuestion();

    return 0;
}