/**
 * File: include/NumCpp.h
 * Author: Ge Yugong
 * Description: 基础数值分析算法库 (No Eigen)
 */

#ifndef NUMCPP_H
#define NUMCPP_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>

// 简单的点结构体
struct Point {
    double t; // Time
    double h; // Height
};

class Interpolator {
public:
    // 拉格朗日插值
    // points: 已知的点集 (例如 t=4,5,9 对应的点)
    // target_t: 需要预测的时间 t
    static double lagrange(const std::vector<Point>& points, double target_t) {
        double result = 0.0;
        int n = points.size();

        for (int i = 0; i < n; ++i) {
            double term = points[i].h; // L_i(x) * y_i
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    // 累乘基函数: (t - t_j) / (t_i - t_j)
                    term *= (target_t - points[j].t) / (points[i].t - points[j].t);
                }
            }
            result += term;
        }
        return result;
    }
};

#endif