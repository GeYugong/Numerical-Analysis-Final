/**
 * File: include/NumCpp.h
 * Author: Ge Yugong
 * Description: 综合数值分析算法库
 * - Part 1: 手写基础算法 (No Eigen, for Problem 1)
 * - Part 2: Eigen 库封装 (For Problem 2)
 */

#ifndef NUMCPP_H
#define NUMCPP_H

// --- 标准库头文件 ---
#include <vector>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <functional>
#include <fstream> // 新增: 用于文件输出
#include <string>  // 新增

// =========================================================
// Part 1: 基础手写算法库 (Problem 1 专用 - 禁用 Eigen)
// =========================================================

// 基础数据结构
struct Point {
    double t; // or x
    double h; // or y
};

// 全局类型定义 (给第一题用)
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

    static double newton(const std::vector<Point>& points, double target_x) {
        int n = points.size();
        Vector divDiff(n);
        for(int i=0; i<n; ++i) divDiff[i] = points[i].h;

        for(int j=1; j<n; ++j) {
            for(int i=n-1; i>=j; --i) {
                divDiff[i] = (divDiff[i] - divDiff[i-1]) / (points[i].t - points[i-j].t);
            }
        }
        double result = divDiff[n-1];
        for(int i=n-2; i>=0; --i) {
            result = result * (target_x - points[i].t) + divDiff[i];
        }
        return result;
    }

    static double piecewiseLinear(const std::vector<Point>& points, double target_x) {
        int n = points.size();
        if (target_x <= points[0].t) return points[0].h;
        if (target_x >= points[n-1].t) return points[n-1].h;
        for(int i=0; i<n-1; ++i) {
            if (target_x >= points[i].t && target_x <= points[i+1].t) {
                double slope = (points[i+1].h - points[i].h) / (points[i+1].t - points[i].t);
                return points[i].h + slope * (target_x - points[i].t);
            }
        }
        return 0.0; 
    }
};

class CubicSpline {
private:
    std::vector<Point> points;
    Vector M; 
public:
    CubicSpline(std::vector<Point> pts) : points(pts) { solveMoments(); }

    void solveMoments() {
        int n = points.size() - 1; 
        if (n < 1) return;
        Vector h(n);
        for(int i=0; i<n; ++i) h[i] = points[i+1].t - points[i].t;
        int size = n - 1;
        if (size == 0) return; 

        Matrix A = LinearAlgebra::zeros(size, size);
        Vector d(size);

        for(int i=0; i<size; ++i) {
            double h_i = h[i];     
            double h_ip1 = h[i+1]; 
            A[i][i] = 2.0 * (h_i + h_ip1);
            if (i > 0) A[i][i-1] = h_i; 
            if (i < size-1) A[i][i+1] = h_ip1; 
            double diff2 = (points[i+2].h - points[i+1].h) / h_ip1;
            double diff1 = (points[i+1].h - points[i].h) / h_i;
            d[i] = 6.0 * (diff2 - diff1);
        }
        Vector M_inner = LinearAlgebra::solveGaussian(A, d);
        M.resize(n+1, 0.0);
        for(int i=0; i<size; ++i) M[i+1] = M_inner[i];
    }

    double evaluate(double x) {
        int n = points.size() - 1;
        int i = 0;
        if (x <= points[0].t) i = 0;
        else if (x >= points[n].t) i = n-1;
        else {
            for(int k=0; k<n; ++k) {
                if (x >= points[k].t && x <= points[k+1].t) {
                    i = k; break;
                }
            }
        }
        double h = points[i+1].t - points[i].t;
        double term1 = M[i] * std::pow(points[i+1].t - x, 3) / (6.0 * h);
        double term2 = M[i+1] * std::pow(x - points[i].t, 3) / (6.0 * h);
        double term3 = (points[i].h - M[i]*h*h/6.0) * (points[i+1].t - x) / h;
        double term4 = (points[i+1].h - M[i+1]*h*h/6.0) * (x - points[i].t) / h;
        return term1 + term2 + term3 + term4;
    }
};

// =========================================================
// Part 2: Eigen 库封装 (Problem 2 专用)
// 使用命名空间 EigenUtils 以避免与 Part 1 的 Vector 定义冲突
// =========================================================
namespace EigenUtils {
    // --- Eigen 头文件 (仅供 Part 2 使用) ---
    // 确保你已经安装了 Eigen 库，并配置了包含路径
    #include <Eigen/Sparse>
    #include <Eigen/Dense>
    #include <Eigen/IterativeLinearSolvers>
    #include <Eigen/SparseCholesky>
    #include <Eigen/SparseLU>
    
    // 类型别名定义
    using SpMat = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using Vector = Eigen::VectorXd; // 注意：这里是 Eigen 的 Vector，不是 std::vector
    using Matrix = Eigen::MatrixXd;

    enum class SolverType { SimplicialLLT, SparseLU };

    class LinearSolver {
    public:
        // 统一求解接口
        static bool solve(const SpMat& A, const Vector& b, Vector& x, SolverType type = SolverType::SimplicialLLT) {
            if (type == SolverType::SimplicialLLT) {
                // 适用于对称正定矩阵 (SPD)，速度最快
                Eigen::SimplicialLLT<SpMat> solver;
                solver.compute(A);
                if(solver.info() != Eigen::Success) {
                    std::cerr << "Decomposition failed!" << std::endl;
                    return false;
                }
                x = solver.solve(b);
                return (solver.info() == Eigen::Success);
            } 
            else if (type == SolverType::SparseLU) {
                // 适用于一般方阵
                Eigen::SparseLU<SpMat> solver;
                solver.compute(A);
                if(solver.info() != Eigen::Success) {
                     std::cerr << "LU Decomposition failed!" << std::endl;
                     return false;
                }
                x = solver.solve(b);
                return (solver.info() == Eigen::Success);
            }
            return false;
        }
    };

    // 导出矩阵到 CSV (用于 Python 可视化)
    inline void saveToCSV(const std::string& filename, const Matrix& matrix) {
        const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream file(filename);
        if (file.is_open()) {
            file << matrix.format(CSVFormat);
            file.close();
            std::cout << "Data saved to " << filename << std::endl;
        } else {
            std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        }
    }
}

#endif