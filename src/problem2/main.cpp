/**
 * File: src/problem2/main.cpp
 * Description: 问题三：均匀四边形网格上的二维稳态热传导 (FVM)
 */

#include <iostream>
#include <vector>
#include "../../include/NumCpp.h" // 确保路径正确

int main() {
    // 1. 物理与网格参数定义
    const int Nx = 50;
    const int Ny = 50;
    const double Lx = 1.0;
    const double Ly = 1.0;
    double dx = Lx / Nx;
    double dy = Ly / Ny;
    double k = 1.0;  // 导热系数
    
    int N_dof = Nx * Ny; // 总自由度 (Degrees of Freedom)

    std::cout << "Initializing FVM Solver..." << std::endl;
    std::cout << "Mesh: " << Nx << "x" << Ny << " (" << N_dof << " cells)" << std::endl;

    // 2. 准备稀疏矩阵的三元组
    std::vector<NumCpp::Triplet> triplets;
    triplets.reserve(N_dof * 5); // 预分配内存，加速填充
    
    NumCpp::Vector b = NumCpp::Vector::Zero(N_dof);

    // 计算传导系数 (Coefficients)
    // Flux = k * Area / Distance
    // East/West Area=dy, Dist=dx -> a_ew = k * dy / dx
    // North/South Area=dx, Dist=dy -> a_ns = k * dx / dy
    double a_e_val = k * dy / dx;
    double a_w_val = k * dy / dx;
    double a_n_val = k * dx / dy;
    double a_s_val = k * dx / dy;

    // 3. 组装全局刚度矩阵与载荷向量
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int P = j * Nx + i; // 当前节点全局索引
            double a_P = 0.0;   // 中心系数 a_P = sum(a_nb)

            // --- West Boundary (i=0) : Dirichlet T=100 ---
            if (i == 0) {
                // 距离减半 (dx/2)，导致通量系数翻倍
                double a_W_bound = 2.0 * a_w_val; 
                a_P += a_W_bound;
                b(P) += a_W_bound * 100.0; // 移项到右端项
            } else {
                // 内部面
                triplets.emplace_back(P, P - 1, -a_w_val);
                a_P += a_w_val;
            }

            // --- East Boundary (i=Nx-1) : Dirichlet T=0 ---
            if (i == Nx - 1) {
                double a_E_bound = 2.0 * a_e_val;
                a_P += a_E_bound;
                b(P) += a_E_bound * 0.0;
            } else {
                triplets.emplace_back(P, P + 1, -a_e_val);
                a_P += a_e_val;
            }

            // --- South Boundary (j=0) : Neumann Adiabatic (Flux=0) ---
            if (j > 0) {
                triplets.emplace_back(P, P - Nx, -a_s_val);
                a_P += a_s_val;
            } 
            // else: j==0, a_S=0, 对 a_P 无贡献 (绝热)

            // --- North Boundary (j=Ny-1) : Neumann Adiabatic (Flux=0) ---
            if (j < Ny - 1) {
                triplets.emplace_back(P, P + Nx, -a_n_val);
                a_P += a_n_val;
            }
            // else: j==Ny-1, a_N=0, 对 a_P 无贡献 (绝热)

            // 填充主对角线
            triplets.emplace_back(P, P, a_P);
        }
    }

    // 4. 构建稀疏矩阵
    NumCpp::SpMat A(N_dof, N_dof);
    A.setFromTriplets(triplets.begin(), triplets.end());

    // 5. 求解线性方程组
    std::cout << "Solving linear system..." << std::endl;
    NumCpp::Vector x(N_dof);
    
    // 使用 SimplicialLLT (Cholesky) 因为矩阵是对称正定的
    bool success = NumCpp::LinearSolver::solve(A, b, x, NumCpp::SolverType::SimplicialLLT);

    if (success) {
        std::cout << "Solver converged successfully!" << std::endl;
        
        // 6. 结果重塑并导出
        // 将一维解向量 x 映射回二维矩阵 (Ny, Nx)
        NumCpp::Matrix T_field(Ny, Nx);
        for(int j=0; j<Ny; ++j) {
            for(int i=0; i<Nx; ++i) {
                T_field(j, i) = x(j * Nx + i);
            }
        }
        
        // 保存到 data 目录
        NumCpp::saveToCSV("data/temperature_field.csv", T_field);
    } else {
        std::cerr << "Solver failed!" << std::endl;
        return 1;
    }

    return 0;
}