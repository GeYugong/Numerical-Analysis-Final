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

int main()
{
    // 1. 读取数据 (注意路径：根据编译执行位置不同，路径可能需要调整)
    // 假设我们在项目根目录执行编译后的程序
    std::string filepath = "D:\\0VS_user\\Numerical-Analysis-Final\\data\\flight_data.csv";
    std::vector<Point> flightData = readCSV(filepath);

    if (flightData.empty())
    {
        std::cerr << "No data loaded. Check path!" << std::endl;
        return 1;
    }

    std::cout << "Data loaded successfully. Total rows: " << flightData.size() << std::endl;

    // 2. 准备插值所需的已知点 (根据题目纠错要求)
    // Group 1: t=4, 5, 9 用于插值 t=6, 7
    std::vector<Point> group1;
    group1.push_back(getPointByTime(flightData, 4.0));
    group1.push_back(getPointByTime(flightData, 5.0));
    group1.push_back(getPointByTime(flightData, 9.0));

    // Group 2: t=5, 9, 10 用于插值 t=8
    std::vector<Point> group2;
    group2.push_back(getPointByTime(flightData, 5.0));
    group2.push_back(getPointByTime(flightData, 9.0));
    group2.push_back(getPointByTime(flightData, 10.0));

    // 3. 执行插值并填补数据
    std::cout << "\n--- Interpolation Results ---" << std::endl;

    // 恢复 t=6
    double h6 = Interpolator::lagrange(group1, 6.0);
    std::cout << "Restored t=6: " << h6 << std::endl;
    // 更新内存中的数据
    flightData[6].h = h6;

    // 恢复 t=7
    double h7 = Interpolator::lagrange(group1, 7.0);
    std::cout << "Restored t=7: " << h7 << std::endl;
    flightData[7].h = h7;

    // 恢复 t=8
    double h8 = Interpolator::lagrange(group2, 8.0);
    std::cout << "Restored t=8: " << h8 << std::endl;
    flightData[8].h = h8;

    // 4. 输出完整数据用于检查
    std::cout << "\n--- Full Data Check ---" << std::endl;
    for (const auto &p : flightData)
    {
        std::cout << "t=" << p.t << ", h=" << p.h << std::endl;
    }

    return 0;
}