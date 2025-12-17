/**
 * File: src/problem1/main.cpp
 * Description: 飞行器数据恢复与分析主程序
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include "../../include/NumCpp.h" // 引用我们在上一层目录创建的头文件

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
            // 简单处理：如果是 NaN，我们标记为 -1 或者在读取逻辑中单独处理
            // 这里我们只存入有效数字，如果读取失败(如NaN字符串)，stod会抛出异常或我们需要手动检查
            // 为了简单，我们只存入我们确定有效的数据行到 data 中，
            // 但考虑到我们需要完整的数组索引，我们这里手动存入全部数据
            // 对于 "NaN"，我们暂时存为 0.0 并在后续逻辑中覆盖
            try
            {
                double h = std::stod(h_str);
                data.push_back({t, h});
            }
            catch (...)
            {
                // 如果转换失败(例如是NaN)，我们存一个标记值
                data.push_back({t, -9999.0});
            }
        }
    }
    return data;
}

// 辅助函数：通过时间t查找高度h (用于构建插值节点)
Point getPointByTime(const std::vector<Point> &all_data, double t)
{
    for (const auto &p : all_data)
    {
        if (std::abs(p.t - t) < 1e-5)
        { // 浮点数比较
            return p;
        }
    }
    return {t, 0.0}; // Should not happen if data is complete
}
// 多项式拟合函数
// 返回系数向量 {c, b, a} 对应 c + bt + at^2
Vector performPolyFit(const std::vector<Point>& data, int degree) {
    int n = data.size();
    int terms = degree + 1; // 二次函数有3个系数

    // 1. 构建设计矩阵 X (Vandermonde matrix)
    // Row i: [1, t_i, t_i^2]
    Matrix X = LinearAlgebra::zeros(n, terms);
    Vector Y(n);

    for (int i = 0; i < n; ++i) {
        double t = data[i].t;
        Y[i] = data[i].h;
        for (int j = 0; j < terms; ++j) {
            X[i][j] = std::pow(t, j);
        }
    }

    // 2. 构建正规方程: (X^T * X) * beta = X^T * Y
    Matrix XT = LinearAlgebra::transpose(X);
    Matrix A = LinearAlgebra::multiply(XT, X);
    Vector B = LinearAlgebra::multiply(XT, Y); // 这里是矩阵乘向量

    // 3. 求解
    Vector coeffs = LinearAlgebra::solveGaussian(A, B);
    return coeffs;
}
int main() {
    // 1. 读取数据 (请保持你的绝对路径)
    std::string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    std::vector<Point> flightData = readCSV(filepath);

    if (flightData.empty()) { std::cerr << "Error loading data.\n"; return 1; }

    // 2. 插值补全 (保留之前的逻辑)
    std::vector<Point> group1 = {getPointByTime(flightData, 4), getPointByTime(flightData, 5), getPointByTime(flightData, 9)};
    std::vector<Point> group2 = {getPointByTime(flightData, 5), getPointByTime(flightData, 9), getPointByTime(flightData, 10)};
    
    flightData[6].h = Interpolator::lagrange(group1, 6.0);
    flightData[7].h = Interpolator::lagrange(group1, 7.0);
    flightData[8].h = Interpolator::lagrange(group2, 8.0);

    std::cout << "Data restoration complete.\n";

    // 3. 执行二次多项式拟合 (h = c + bt + at^2)
    // 我们使用所有数据点进行拟合
    std::cout << "\n--- Performing Quadratic Regression ---" << std::endl;
    Vector coeffs = performPolyFit(flightData, 2);

    double c = coeffs[0];
    double b = coeffs[1];
    double a = coeffs[2];

    std::cout << "Fitted Model: h(t) = " << a << " * t^2 + " << b << " * t + " << c << std::endl;

    // 4. 计算均方误差 (MSE)
    double mse = 0.0;
    for (const auto& p : flightData) {
        double pred = a * p.t * p.t + b * p.t + c;
        mse += std::pow(pred - p.h, 2);
    }
    mse /= flightData.size();
    std::cout << "MSE: " << mse << std::endl;

    // 5. 导出拟合结果到CSV (用于画图)
    // 我们可以生成更密集的点来画平滑曲线
    std::ofstream outfile("data/regression_result.csv");
    outfile << "t,h_real,h_pred\n";
    for(const auto& p : flightData) {
        double pred = a * p.t * p.t + b * p.t + c;
        outfile << p.t << "," << p.h << "," << pred << "\n";
    }
    std::cout << "Regression data saved to data/regression_result.csv" << std::endl;

    return 0;
}