/**
 * File: include/NumCpp.h
 * Author: Ge Yugong
 * Description: 基础数值分析算法库 (No Eigen) - 包含插值(Lagrange/Newton/Spline)、矩阵求解、非线性方程
 */

#ifndef NUMCPP_H
#define NUMCPP_H

#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <functional>

// 基础数据结构
struct Point {
    double t; // or x
    double h; // or y
};

using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

class LinearAlgebra {
public:
    static Matrix zeros(int rows, int cols) { return Matrix(rows, Vector(cols, 0.0)); }

    static Matrix transpose(const Matrix& A) {
        if (A.empty()) return {};
        int rows = A.size(), cols = A[0].size();
        Matrix res = zeros(cols, rows);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j) res[j][i] = A[i][j];
        return res;
    }

    static Matrix multiply(const Matrix& A, const Matrix& B) {
        if (A.empty() || B.empty() || A[0].size() != B.size())
            throw std::invalid_argument("Matrix dimension mismatch");
        int rows = A.size(), cols = B[0].size(), common = B.size();
        Matrix res = zeros(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                for (int k = 0; k < common; ++k) res[i][j] += A[i][k] * B[k][j];
        return res;
    }

    static Vector multiply(const Matrix& A, const Vector& x) {
        if (A.empty() || A[0].size() != x.size())
            throw std::invalid_argument("Dimension mismatch");
        int rows = A.size(), cols = x.size();
        Vector res(rows, 0.0);
        for(int i=0; i<rows; ++i)
            for(int j=0; j<cols; ++j) res[i] += A[i][j] * x[j];
        return res;
    }

    static Vector solveGaussian(Matrix A, Vector b) {
        int n = A.size();
        for (int i = 0; i < n; ++i) {
            int pivot = i;
            for (int k = i + 1; k < n; ++k)
                if (std::abs(A[k][i]) > std::abs(A[pivot][i])) pivot = k;
            std::swap(A[i], A[pivot]);
            std::swap(b[i], b[pivot]);
            
            // Singular check omitted for brevity in robust implementation
            
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; ++j) A[k][j] -= factor * A[i][j];
                b[k] -= factor * b[i];
            }
        }
        Vector x(n);
        for (int i = n - 1; i >= 0; --i) {
            double sum = 0.0;
            for (int j = i + 1; j < n; ++j) sum += A[i][j] * x[j];
            x[i] = (b[i] - sum) / A[i][i];
        }
        return x;
    }

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
                for (int j = 0; j < n; ++j) if (j != i) sigma += A[i][j] * x_new[j];
                double x_sor = (1 - omega) * x[i] + (omega / A[i][i]) * (b[i] - sigma);
                if (std::abs(x_sor - x[i]) > max_diff) max_diff = std::abs(x_sor - x[i]);
                x_new[i] = x_sor;
            }
            x = x_new;
            if (max_diff < tol) break;
        }
        return {x, iter + 1};
    }
};

class Solvers {
public:
    static double newtonMethod(std::function<double(double)> func, 
                               std::function<double(double)> deriv, 
                               double x0, double tol = 1e-6, int max_iter = 100) {
        double x = x0;
        for (int i = 0; i < max_iter; ++i) {
            double f_val = func(x);
            double f_deriv = deriv(x);
            if (std::abs(f_deriv) < 1e-10) break; 
            double step = f_val / f_deriv;
            x = x - step;
            if (std::abs(step) < tol) return x;
        }
        return x;
    }
};

class Interpolator {
public:
    // 1. Lagrange Interpolation
    static double lagrange(const std::vector<Point>& points, double target_x) {
        double result = 0.0;
        int n = points.size();
        for (int i = 0; i < n; ++i) {
            double term = points[i].h;
            for (int j = 0; j < n; ++j)
                if (i != j) term *= (target_x - points[j].t) / (points[i].t - points[j].t);
            result += term;
        }
        return result;
    }

    // 2. Newton Interpolation (Divided Differences)
    static double newton(const std::vector<Point>& points, double target_x) {
        int n = points.size();
        // 计算均差表 (Divided Difference Table)
        // 为了节省空间，我们只需要一维数组存储对角线元素
        Vector divDiff(n);
        for(int i=0; i<n; ++i) divDiff[i] = points[i].h;

        for(int j=1; j<n; ++j) {
            for(int i=n-1; i>=j; --i) {
                divDiff[i] = (divDiff[i] - divDiff[i-1]) / (points[i].t - points[i-j].t);
            }
        }

        // 计算多项式值: N(x) = c0 + c1(x-x0) + c2(x-x0)(x-x1) + ...
        double result = divDiff[n-1];
        for(int i=n-2; i>=0; --i) {
            result = result * (target_x - points[i].t) + divDiff[i];
        }
        return result;
    }

    // 3. Piecewise Linear Interpolation
    static double piecewiseLinear(const std::vector<Point>& points, double target_x) {
        // 假设 points 已按 t 排序
        int n = points.size();
        if (target_x <= points[0].t) return points[0].h;
        if (target_x >= points[n-1].t) return points[n-1].h;

        // 寻找区间 [t_i, t_{i+1}]
        for(int i=0; i<n-1; ++i) {
            if (target_x >= points[i].t && target_x <= points[i+1].t) {
                double slope = (points[i+1].h - points[i].h) / (points[i+1].t - points[i].t);
                return points[i].h + slope * (target_x - points[i].t);
            }
        }
        return 0.0; 
    }
};

// --- 新增: 三次样条插值类 (需要存储状态) ---
class CubicSpline {
private:
    std::vector<Point> points;
    Vector M; // 二阶导数值

public:
    CubicSpline(std::vector<Point> pts) : points(pts) {
        solveMoments();
    }

    // 构建三对角矩阵并求解 M (自然边界条件 M0=Mn=0)
    void solveMoments() {
        int n = points.size() - 1; // n个区间，n+1个点
        if (n < 1) return;

        Vector h(n);
        for(int i=0; i<n; ++i) h[i] = points[i+1].t - points[i].t;

        // 构建线性方程组 A * M = d
        // 矩阵大小 (n-1) * (n-1) -> 对应 M_1 到 M_{n-1}
        int size = n - 1;
        if (size == 0) return; // 只有2个点，退化为线性，M均为0

        Matrix A = LinearAlgebra::zeros(size, size);
        Vector d(size);

        for(int i=0; i<size; ++i) {
            // 对应内部节点 i+1
            double h_i = h[i];     // h_{i}
            double h_ip1 = h[i+1]; // h_{i+1}
            
            // 对角线元素
            A[i][i] = 2.0 * (h_i + h_ip1);
            
            // 邻居
            if (i > 0) A[i][i-1] = h_i; // mu
            if (i < size-1) A[i][i+1] = h_ip1; // lambda

            // 右端项 d_i (利用均差)
            double diff2 = (points[i+2].h - points[i+1].h) / h_ip1;
            double diff1 = (points[i+1].h - points[i].h) / h_i;
            d[i] = 6.0 * (diff2 - diff1);
        }

        Vector M_inner = LinearAlgebra::solveGaussian(A, d);

        // 填充完整的 M (包含边界 0)
        M.resize(n+1, 0.0);
        for(int i=0; i<size; ++i) M[i+1] = M_inner[i];
    }

    double evaluate(double x) {
        int n = points.size() - 1;
        // 找到 x 所在的区间 [x_i, x_{i+1}]
        int i = 0;
        if (x <= points[0].t) i = 0;
        else if (x >= points[n].t) i = n-1;
        else {
            for(int k=0; k<n; ++k) {
                if (x >= points[k].t && x <= points[k+1].t) {
                    i = k;
                    break;
                }
            }
        }

        double h = points[i+1].t - points[i].t;
        double a = points[i].t;
        // 样条公式 S(x)
        double term1 = M[i] * std::pow(points[i+1].t - x, 3) / (6.0 * h);
        double term2 = M[i+1] * std::pow(x - points[i].t, 3) / (6.0 * h);
        double term3 = (points[i].h - M[i]*h*h/6.0) * (points[i+1].t - x) / h;
        double term4 = (points[i+1].h - M[i+1]*h*h/6.0) * (x - points[i].t) / h;

        return term1 + term2 + term3 + term4;
    }
};

#endif