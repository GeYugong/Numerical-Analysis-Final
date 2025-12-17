/**
 * @file NumUtils.h
 * @brief 基于 Eigen 库的数值分析通用工具集
 */
#ifndef NUM_UTILS_H
#define NUM_UTILS_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

namespace NumUtils { // 命名空间已重命名
    using SpMat = Eigen::SparseMatrix<double>;
    using Triplet = Eigen::Triplet<double>;
    using Vector = Eigen::VectorXd;
    using Matrix = Eigen::MatrixXd;

    enum class SolverType { SimplicialLLT, SparseLU };

    class LinearSolver {
    public:
        // 求解 Ax = b
        static bool solve(const SpMat& A, const Vector& b, Vector& x, SolverType type = SolverType::SimplicialLLT) {
            if (type == SolverType::SimplicialLLT) {
                // SimplicialLLT 适用于对称正定矩阵 (SPD)，本题的热传导矩阵符合此特性
                Eigen::SimplicialLLT<SpMat> solver(A);
                x = solver.solve(b);
                return (solver.info() == Eigen::Success);
            } 
            else if (type == SolverType::SparseLU) {
                // SparseLU 适用于一般方阵
                Eigen::SparseLU<SpMat> solver;
                solver.compute(A);
                x = solver.solve(b);
                return (solver.info() == Eigen::Success);
            }
            return false;
        }
    };

    // 将矩阵数据导出为 CSV，方便 Python 读取
    inline void saveToCSV(const std::string& filename, const Matrix& matrix) {
        const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");
        std::ofstream file(filename);
        if (file.is_open()) {
            file << matrix.format(CSVFormat);
            file.close();
            std::cout << "Data successfully saved to " << filename << std::endl;
        } else {
            std::cerr << "Error: Unable to open file " << filename << std::endl;
        }
    }
}
#endif