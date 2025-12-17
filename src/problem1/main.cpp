/**
 * File: src/problem1/main.cpp
 * Description: 飞行器数据恢复与分析主程序 (Final Version)
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

// 辅助函数：读取CSV文件
std::vector<Point> readCSV(const std::string &filename)
{
    std::vector<Point> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
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
        if (std::abs(p.t - t) < 1e-5) return p;
    }
    return {t, 0.0}; 
}

// 多项式拟合函数
Vector performPolyFit(const std::vector<Point>& data, int degree) {
    int n = data.size();
    int terms = degree + 1; 

    // 1. 构建设计矩阵 X
    Matrix X = LinearAlgebra::zeros(n, terms);
    Vector Y(n);

    for (int i = 0; i < n; ++i) {
        double t = data[i].t;
        Y[i] = data[i].h;
        for (int j = 0; j < terms; ++j) {
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

int main() {
    // ==========================================
    // Part 1: 数据读取与插值恢复
    // ==========================================
    std::string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    std::vector<Point> flightData = readCSV(filepath);

    if (flightData.empty()) { std::cerr << "Error loading data.\n"; return 1; }

    std::vector<Point> group1 = {getPointByTime(flightData, 4), getPointByTime(flightData, 5), getPointByTime(flightData, 9)};
    std::vector<Point> group2 = {getPointByTime(flightData, 5), getPointByTime(flightData, 9), getPointByTime(flightData, 10)};
    
    // 补全数据
    flightData[6].h = Interpolator::lagrange(group1, 6.0);
    flightData[7].h = Interpolator::lagrange(group1, 7.0);
    flightData[8].h = Interpolator::lagrange(group2, 8.0);

    std::cout << "Data restoration complete.\n";

    // ==========================================
    // Part 2: 二次多项式回归
    // ==========================================
    std::cout << "\n--- [Step 2] Performing Quadratic Regression ---" << std::endl;
    Vector coeffs = performPolyFit(flightData, 2); // 结果为 [c, b, a]

    double c = coeffs[0];
    double b = coeffs[1];
    double a = coeffs[2];

    std::cout << "Fitted Model: h(t) = " << a << " * t^2 + " << b << " * t + " << c << std::endl;

    // 计算MSE
    double mse = 0.0;
    for (const auto& p : flightData) {
        double pred = a * p.t * p.t + b * p.t + c;
        mse += std::pow(pred - p.h, 2);
    }
    mse /= flightData.size();
    std::cout << "MSE: " << mse << std::endl;

    // 导出拟合结果
    std::ofstream outfile("data/regression_result.csv");
    outfile << "t,h_real,h_pred\n";
    for(const auto& p : flightData) {
        double pred = a * p.t * p.t + b * p.t + c;
        outfile << p.t << "," << p.h << "," << pred << "\n";
    }
    std::cout << "Regression data saved to data/regression_result.csv" << std::endl;

    // ==========================================
    // Part 3: 落地时间预测 (Newton法)
    // ==========================================
    std::cout << "\n--- [Step 3] Landing Time Prediction ---" << std::endl;
    
    // 定义高度函数 h(t) 和其导数 h'(t) = 2at + b
    auto h_func = [&](double t) { return a * t * t + b * t + c; };
    auto h_deriv = [&](double t) { return 2 * a * t + b; };

    // 使用牛顿法，初始猜测 t=14
    double landing_time = Solvers::newtonMethod(h_func, h_deriv, 14.0);
    
    std::cout << "Predicted Landing Time (h=0): " << std::setprecision(6) << std::fixed << landing_time << " s" << std::endl;

    // ==========================================
    // Part 4: 引擎推力解算 (SOR迭代法)
    // ==========================================
    std::cout << "\n--- [Step 4] Engine Thrust Control (SOR Method) ---" << std::endl;
    
    // 方程组:
    // 4x1 + 3x2       = 24
    // 3x1 + 4x2 - x3  = 30
    //       -x2 + 4x3 = -24
    
    Matrix A_eng = {
        {4.0,  3.0,  0.0},
        {3.0,  4.0, -1.0},
        {0.0, -1.0,  4.0}
    };
    Vector b_eng = {24.0, 30.0, -24.0};
    
    Vector x0 = {1.0, 1.0, 1.0}; // 初始猜测
    double omega = 1.25;         // 松弛因子
    double tol = 1e-6;           // 容差
    
    // 调用 SOR 求解
    std::pair<Vector, int> result = LinearAlgebra::solveSOR(A_eng, b_eng, omega, x0, tol, 100);
    Vector thrusts = result.first;
    int iterations = result.second;

    std::cout << "SOR Converged in " << iterations << " iterations." << std::endl;
    std::cout << "Engine Thrusts Solution:" << std::endl;
    std::cout << "x1 = " << thrusts[0] << std::endl;
    std::cout << "x2 = " << thrusts[1] << std::endl;
    std::cout << "x3 = " << thrusts[2] << std::endl;

    return 0;
}