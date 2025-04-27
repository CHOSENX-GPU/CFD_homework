/////////////////////////////////////////////////////////////////////////////
// Module.h - 全局变量与参数定义
/////////////////////////////////////////////////////////////////////////////

#pragma once
#include <vector>
#include <memory>
#include <string>

using namespace std;

// 模拟参数命名空间
namespace CFDParams
{
    // 计算网格参数
    extern int ni;                 // 网格节点数
    extern double x1;              // 起始坐标
    extern double x2;              // 终止坐标
    extern double dx;              // 网格步长
    extern unique_ptr<double[]> x; // 网格坐标

    // 计算场变量
    extern unique_ptr<double[]> u;  // 当前时刻速度
    extern unique_ptr<double[]> um; // n+1时刻速度
    extern unique_ptr<double[]> un; // n-1时刻速度

    // 时间参数
    extern double t0;   // 当前时间
    extern double tout; // 终止时间
    extern double dt;   // 时间步长

    // 物理参数
    extern double cfl; // CFL数
    extern double a;   // 波速

    // 输出控制
    extern int nloop; // 当前迭代步数
    extern int nout;  // 输出间隔

    // 求解方法枚举
    enum SolverMethod
    {
        UPWIND = 1,            // 迎风格式
        FTCS = 2,              // 时间向前空间中心差分
        LAX_FRIEDRICHS = 3,    // Lax-Friedrichs格式
        LAX_WENDROFF = 4,      // Lax-Wendroff格式
        LEAPFROG = 5,          // Leap-frog格式
        MACCORMACK = 6,        // MacCormack格式
        ROE = 7,               // Roe格式
        BEAM_WARMING = 8,      // Beam-Warming显式格式
        IMPLICIT_UPWIND = 9,   // 隐式迎风格式
        CRANK_NICOLSON = 10,   // Crank-Nicolson格式
        BEAM_WARMING_IMP = 11, // Beam-Warming隐式格式
        RK3_WENO5 = 12         // 三阶Runge-Kutta + 五阶WENO格式
    };

    extern SolverMethod currentMethod; // 当前使用的求解方法

    // 获取求解方法名称
    string getSolverName(SolverMethod method);

    // 初始化全局参数
    bool initialize();
}