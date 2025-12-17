/**
 * File: include/NumCpp.h
 * Author: Ge Yugong
 * Description: 基础数值分析算法库 (No Eigen) - 包含插值、矩阵求解(高斯/SOR)、非线性方程求解
 */

#ifndef NUMCPP_H
#define NUMCPP_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <functional> // 新增：用于传递数学函数

// 基础数据结构
struct Point {
    double t;
    double h;
};

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

class LinearAlgebra {
public:
    // 创建一个 rows x cols 的零矩阵
    static Matrix zeros(int rows, int cols) {
        return Matrix(rows, Vector(cols, 0.0));
    }

    // 矩阵转置 (A^T)
    static Matrix transpose(const Matrix& A) {
        if (A.empty()) return {};
        int rows = A.size();
        int cols = A[0].size();
        Matrix res = zeros(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                res[j][i] = A[i][j];
        return res;
    }

    // 矩阵乘法 (A * B)
    static Matrix multiply(const Matrix& A, const Matrix& B) {
        if (A.empty() || B.empty() || A[0].size() != B.size())
            throw std::invalid_argument("Matrix dimension mismatch for multiplication");
        
        int rows = A.size();
        int cols = B[0].size();
        int common = B.size();
        Matrix res = zeros(rows, cols);

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                for (int k = 0; k < common; ++k) {
                    res[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return res;
    }

    // 矩阵乘以向量 (A * x)
    static Vector multiply(const Matrix& A, const Vector& x) {
        if (A.empty() || A[0].size() != x.size())
            throw std::invalid_argument("Dimension mismatch for Mat-Vec mul");
        
        int rows = A.size();
        int cols = x.size();
        Vector res(rows, 0.0);
        
        for(int i=0; i<rows; ++i) {
            for(int j=0; j<cols; ++j) {
                res[i] += A[i][j] * x[j];
            }
        }
        return res;
    }

    // 高斯消元法求解线性方程组 Ax = b
    static Vector solveGaussian(Matrix A, Vector b) {
        int n = A.size();
        if (n == 0 || A[0].size() != n || b.size() != n)
            throw std::invalid_argument("Invalid system for Gaussian elimination");

        // 1. 前向消元
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[pivot][i])) pivot = k;
            }
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);

            if (std::abs(A[i][i]) < 1e-10) 
                throw std::runtime_error("Matrix is singular!");

            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }

        // 2. 回代
        Vector x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

    // --- 新增功能: SOR 迭代求解器 ---
    // 返回 pair<解向量, 迭代次数>
    static std::pair<Vector, int> solveSOR(const Matrix& A, const Vector& b, 
                                           double omega, Vector x0, 
                                           double tol, int max_iter) {
        int n = A.size();
        Vector x = x0;
        int iter = 0;
        
        for (; iter < max_iter; ++iter) {
            Vector x_new = x;
            double max_diff = 0.0;
            
            for (int i = 0; i < n; ++i) {
                double sigma = 0.0;
                // 利用 Gauss-Seidel 特性: j < i 使用新值, j > i 使用旧值
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        // 注意：这里 x_new[j] 在 j<i 时已经是更新过的值了
                        sigma += A[i][j] * x_new[j]; 
                    }
                }
                // SOR 公式: x_new = (1-w)x_old + (w/aii)(b - sigma)
                double x_sor = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma);
                
                double diff = std::abs(x_sor - x[i]);
                if (diff > max_diff) max_diff = diff;
                
                x_new[i] = x_sor;
            }
            x = x_new;
            if (max_diff < tol) break;
        }
        return {x, iter + 1};
    }
};

// --- 新增类: 非线性方程求解器 ---
class Solvers {
public:
    // 牛顿迭代法求根: x_{k+1} = x_k - f(x)/f'(x)
    static double newtonMethod(std::function<double(double)> func, 
                               std::function<double(double)> deriv, 
                               double x0, double tol = 1e-6, int max_iter = 100) {
        double x = x0;
        for (int i = 0; i < max_iter; ++i) {
            double f_val = func(x);
            double f_deriv = deriv(x);
            
            if (std::abs(f_deriv) < 1e-10) {
                std::cerr << "Warning: Derivative too close to 0 in Newton method.\n";
                break; 
            }
            
            double step = f_val / f_deriv;
            x = x - step;
            
            if (std::abs(step) < tol) return x;
        }
        return x;
    }
};

class Interpolator {
public:
    static double lagrange(const std::vector<Point>& points, double target_t) {
        double result = 0.0;
        int n = points.size();
        for (int i = 0; i < n; ++i) {
            double term = points[i].h;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    term *= (target_t - points[j].t) / (points[i].t - points[j].t);
                }
            }
            result += term;
        }
        return result;
    }
};

#endif