///////////////////////////////////////////////////////////////////
// Solver.cpp - 数值求解算法库
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <algorithm>
#include "Module.h"

using namespace CFDParams;
using namespace std;

// 辅助函数：获取指定符号
inline double sign(double x) {
    return (x > 0.0) ? 1.0 : -1.0;
}

// 边界点处理函数
inline double bound(int i, const double* u_array) {
    if (i < 0) return u_array[0];     // 左边界外推
    if (i >= ni) return u_array[ni-1]; // 右边界外推
    return u_array[i];                // 内部点
}

// 五阶WENO重构
double WENO5(double vmm, double vm, double v, double vp, double vpp) {
    const double eps = 1.0e-6; // 避免除以零
    
    // 平滑度指标
    double S1 = 13.0/12.0 * pow(vmm - 2.0*vm + v, 2) + 
                1.0/4.0 * pow(vmm - 4.0*vm + 3.0*v, 2);
    
    double S2 = 13.0/12.0 * pow(vm - 2.0*v + vp, 2) + 
                1.0/4.0 * pow(vm - vp, 2);
    
    double S3 = 13.0/12.0 * pow(v - 2.0*vp + vpp, 2) + 
                1.0/4.0 * pow(3.0*v - 4.0*vp + vpp, 2);
    
    // 各模板的理想权重
    const double gamma1 = 0.1;
    const double gamma2 = 0.6;
    const double gamma3 = 0.3;
    
    // 非归一化权重
    double alpha1 = gamma1 / pow(eps + S1, 2);
    double alpha2 = gamma2 / pow(eps + S2, 2);
    double alpha3 = gamma3 / pow(eps + S3, 2);
    
    // 归一化权重
    double alpha_sum = alpha1 + alpha2 + alpha3;
    double w1 = alpha1 / alpha_sum;
    double w2 = alpha2 / alpha_sum;
    double w3 = alpha3 / alpha_sum;
    
    // 各模板的多项式近似
    double p1 = (1.0/3.0)*vmm - (7.0/6.0)*vm + (11.0/6.0)*v;
    double p2 = -(1.0/6.0)*vm + (5.0/6.0)*v + (1.0/3.0)*vp;
    double p3 = (1.0/3.0)*v + (5.0/6.0)*vp - (1.0/6.0)*vpp;
    
    // 最终重构结果
    return w1*p1 + w2*p2 + w3*p3;
}

// 计算数值通量
double computeFlux(int i, const double* u_array) {
    double fmm, fm, f, fp, fpp;
    
    // 根据波速方向选择上风格式
    if (a < 0) { // 波速为负，右侧为上风向
        fmm = (bound(i+3, u_array) - bound(i+2, u_array)) / dx;
        fm = (bound(i+2, u_array) - bound(i+1, u_array)) / dx;
        f = (bound(i+1, u_array) - bound(i, u_array)) / dx;
        fp = (bound(i, u_array) - bound(i-1, u_array)) / dx;
        fpp = (bound(i-1, u_array) - bound(i-2, u_array)) / dx;
    } else { // 波速为正，左侧为上风向
        fmm = (bound(i-2, u_array) - bound(i-3, u_array)) / dx;
        fm = (bound(i-1, u_array) - bound(i-2, u_array)) / dx;
        f = (bound(i, u_array) - bound(i-1, u_array)) / dx;
        fp = (bound(i+1, u_array) - bound(i, u_array)) / dx;
        fpp = (bound(i+2, u_array) - bound(i+1, u_array)) / dx;
    }
    
    // 使用WENO5重构计算导数
    return -a * WENO5(fmm, fm, f, fp, fpp);
}

// 边界条件处理
void call_boundary() {
    u[0] = 0.0;    // 左边界速度为0
    u[ni-1] = 0.0; // 右边界速度为0
}

// 基于CFL条件计算时间步长
void call_CFL() {
    dt = cfl * dx / abs(a);
    cout << "时间步长dt = " << dt << endl;
}

// 格式1: 一阶迎风格式
void call_solve_1() {
    for (int i = 1; i < ni-1; i++) {
        if (a > 0) {
            // 正波速，使用向后差分
            um[i] = u[i] - a * dt / dx * (u[i] - u[i-1]);
        } else {
            // 负波速，使用向前差分
            um[i] = u[i] - a * dt / dx * (u[i+1] - u[i]);
        }
    }
    
    // 更新解
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式2: 时间向前空间中心差分格式（FTCS）
void call_solve_2() {
    for (int i = 1; i < ni-1; i++) {
        um[i] = u[i] - a * dt / (2 * dx) * (u[i+1] - u[i-1]);
    }

    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式3: Lax-Friedrichs格式
void call_solve_3() {
    for (int i = 1; i < ni-1; i++) {
        um[i] = 0.5 * (u[i+1] + u[i-1]) - 0.5 * a * dt / dx * (u[i+1] - u[i-1]);
    }

    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式4: Lax-Wendroff格式
void call_solve_4() {
    for (int i = 1; i < ni-1; i++) {
        um[i] = u[i] - 0.5 * a * dt / dx * (u[i+1] - u[i-1]) + 
                0.5 * a * a * dt * dt / (dx * dx) * (u[i+1] - 2 * u[i] + u[i-1]);
    }

    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式5: Leap-frog格式
void call_solve_5() {
    if (nloop == 0) {
        // 第一步用迎风格式初始化
        for (int i = 1; i < ni-1; i++) {
            um[i] = u[i] - a * dt / dx * (u[i] - u[i-1]);
        }
    } else {
        for (int i = 1; i < ni-1; i++) {
            um[i] = un[i] - a * dt / dx * (u[i+1] - u[i-1]);
        }
    }
    
    for (int i = 1; i < ni-1; i++) {
        un[i] = u[i];
        u[i] = um[i];
    }
}

// 格式6: MacCormack格式
void call_solve_6() {
    // 预测步
    for (int i = 1; i < ni-1; i++) {
        um[i] = u[i] - a * dt / dx * (u[i+1] - u[i]);
    }
    
    // 校正步
    for (int i = 1; i < ni-1; i++) {
        u[i] = 0.5 * (u[i] + um[i]) - 0.5 * a * dt / dx * (um[i] - um[i-1]);
    }
}

// 格式8: Beam-Warming显式格式
void call_solve_7() {
    um[1] = u[1] - a * dt / dx * (u[1] - u[0]); // 处理第二个点
    
    for (int i = 2; i < ni-1; i++) {
        um[i] = u[i] - 0.5 * a * dt / dx * (3 * u[i] - 4 * u[i-1] + u[i-2]);
    }
    
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式9: 隐式迎风格式
void call_solve_8() {
    double a1 = a * dt / dx;
    double a2 = 1.0 + a1;
    
    um[0] = u[0]; // 处理边界
    
    for (int i = 1; i < ni-1; i++) {
        um[i] = (u[i] + a1 * um[i-1]) / a2;
    }
    
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式10: Crank-Nicolson格式
void call_solve_9() {
    double a1 = a * dt / (2 * dx);
    double a2 = 1.0 + a1;
    
    um[0] = u[0]; // 处理边界
    
    for (int i = 1; i < ni-1; i++) {
        um[i] = (u[i] + a1 * (um[i-1] + u[i-1] - u[i])) / a2;
    }
    
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式11: Beam-Warming隐式格式
void call_solve_10() {
    um[0] = u[0]; // 处理边界
    
    // 处理第二个点
    um[1] = (u[1] + a * dt / (2 * dx) * (um[0] + u[0] - u[1])) / (1.0 + a * dt / (2 * dx));
    
    double a1 = -a * dt / (2 * dx);
    double a2 = 2 * a * dt / dx;
    double a3 = 1.0 + a1 + a2;
    
    for (int i = 2; i < ni-1; i++) {
        um[i] = (u[i] + a1 * um[i-2] + a2 * um[i-1]) / a3;
    }
    
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 格式12: 三阶Runge-Kutta + 五阶WENO格式
void call_solve_11() {
    // 分配临时数组
    auto u1 = std::make_unique<double[]>(ni);
    auto u2 = std::make_unique<double[]>(ni);
    
    // 第一步 RK3
    for (int i = 1; i < ni-1; i++) {
        u1[i] = u[i] + dt * computeFlux(i, u.get());
    }
    u1[0] = u[0];        // 边界
    u1[ni-1] = u[ni-1];  // 边界
    
    // 第二步 RK3
    for (int i = 1; i < ni-1; i++) {
        u2[i] = 0.75 * u[i] + 0.25 * u1[i] + 0.25 * dt * computeFlux(i, u1.get());
    }
    u2[0] = u[0];        // 边界
    u2[ni-1] = u[ni-1];  // 边界
    
    // 第三步 RK3
    for (int i = 1; i < ni-1; i++) {
        um[i] = (1.0/3.0) * u[i] + (2.0/3.0) * u2[i] + 
                (2.0/3.0) * dt * computeFlux(i, u2.get());
    }
    
    // 更新解
    for (int i = 1; i < ni-1; i++) {
        u[i] = um[i];
    }
}

// 求解器选择函数
void callSolver(SolverMethod method) {
    switch (method) {
        case UPWIND:
            call_solve_1();
            break;
        case FTCS:
            call_solve_2();
            break;
        case LAX_FRIEDRICHS:
            call_solve_3();
            break;
        case LAX_WENDROFF:
            call_solve_4();
            break;
        case LEAPFROG:
            call_solve_5();
            break;
        case MACCORMACK:
            call_solve_6();
            break;
        case BEAM_WARMING:
            call_solve_7();
            break;
        case IMPLICIT_UPWIND:
            call_solve_8();
            break;
        case CRANK_NICOLSON:
            call_solve_9();
            break;
        case BEAM_WARMING_IMP:
            call_solve_10();
            break;
        case RK3_WENO5:
            call_solve_11();
            break;
        default:
            cout << "未知求解方法，使用默认的RK3_WENO5格式" << endl;
            call_solve_11();
            break;
    }
}