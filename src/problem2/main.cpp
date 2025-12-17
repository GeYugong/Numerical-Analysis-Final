#include "../../include/NumUtils.h"
#include <iostream>

int main() {
    // 1. 参数定义 (对应题目要求：50*50均匀网格 [cite: 36])
    const int Nx = 50, Ny = 50;
    double dx = 1.0/Nx, dy = 1.0/Ny, k = 1.0;
    int N_dof = Nx * Ny;

    std::cout << "Initializing FVM Solver with Grid: " << Nx << "x" << Ny << "..." << std::endl;

    // 2. 准备稀疏矩阵的三元组
    std::vector<NumUtils::Triplet> triplets;
    // 预分配内存以提升填充速度 (5点差分格式，每行最多5个非零元)
    triplets.reserve(N_dof * 5); 
    
    NumUtils::Vector b = NumUtils::Vector::Zero(N_dof);

    // 系数计算 (对应离散方程系数)
    double a_e = k * dy / dx;
    double a_w = k * dy / dx;
    double a_n = k * dx / dy;
    double a_s = k * dx / dy;

    // 3. 组装全局刚度矩阵与负载向量
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int P = j * Nx + i; // 当前节点索引
            double a_P = 0.0;   // 中心系数累加器

            // --- 处理西边界 (Left, i=0) T=100 ---
            if (i == 0) { 
                // Dirichlet 边界修正: a_w 项移至右端项 b
                double dist_correction = 2.0; // 距离为 h/2，通量增加一倍
                a_P += dist_correction * a_w; 
                b(P) += dist_correction * a_w * 100.0; 
            } else { 
                triplets.emplace_back(P, P-1, -a_w); 
                a_P += a_w; 
            }

            // --- 处理东边界 (Right, i=Nx-1) T=0 ---
            if (i == Nx-1) { 
                // Dirichlet 边界修正: T=0，右端项增加 0
                double dist_correction = 2.0;
                a_P += dist_correction * a_e; 
                b(P) += dist_correction * a_e * 0.0; 
            } else { 
                triplets.emplace_back(P, P+1, -a_e); 
                a_P += a_e; 
            }

            // --- 处理南边界 (Bottom, j=0) 绝热 ---
            if (j > 0) { 
                triplets.emplace_back(P, P-Nx, -a_s); 
                a_P += a_s; 
            }
            // 绝热边界处通量为0，不需要对 a_P 或 b 做额外操作，相当于 a_s=0

            // --- 处理北边界 (Top, j=Ny-1) 绝热 ---
            if (j < Ny-1) { 
                triplets.emplace_back(P, P+Nx, -a_n); 
                a_P += a_n; 
            }

            // 对角线主元素
            triplets.emplace_back(P, P, a_P);
        }
    }

    // 4. 构建稀疏矩阵并求解
    NumUtils::SpMat A(N_dof, N_dof);
    A.setFromTriplets(triplets.begin(), triplets.end());
    
    NumUtils::Vector x(N_dof);
    std::cout << "Solving linear system..." << std::endl;
    
    // 使用 LLT 分解 (Cholesky)，速度最快
    if(NumUtils::LinearSolver::solve(A, b, x, NumUtils::SolverType::SimplicialLLT)) {
        std::cout << "Solved successfully!" << std::endl;
    } else {
        std::cerr << "Solver failed!" << std::endl;
        return -1;
    }

    // 5. 结果重构与导出
    NumUtils::Matrix T_field(Ny, Nx);
    for(int j=0; j<Ny; ++j) {
        for(int i=0; i<Nx; ++i) {
            // 注意：存储时通常行优先，对应物理坐标 y, x
            T_field(j,i) = x(j*Nx+i);
        }
    }
    
    NumUtils::saveToCSV("temperature_field.csv", T_field);
    return 0;
}