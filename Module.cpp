//////////////////////////////////////////////////////////////////////////
// Module.cpp - 全局变量定义与初始化
//////////////////////////////////////////////////////////////////////////

#include "Module.h"
#include <iostream>
using namespace std;
namespace CFDParams
{
    // 网格参数
    int ni = 101;
    double x1 = 0.0;
    double x2 = 1.0;
    double dx = 0.0;

    // 计算变量
    unique_ptr<double[]> x = nullptr;
    unique_ptr<double[]> u = nullptr;
    unique_ptr<double[]> um = nullptr;
    unique_ptr<double[]> un = nullptr;

    // 时间参数
    double t0 = 0.0;
    double tout = 3.0;
    double dt = 0.0;

    // 物理参数
    double cfl = 0.2;
    double a = 0.2;

    // 输出控制
    int nloop = 0;
    int nout = 2;

    // 默认求解方法
    SolverMethod currentMethod = RK3_WENO5;

    // 获取求解方法名称
    string getSolverName(SolverMethod method)
    {
        switch (method)
        {
        case UPWIND:
            return "Upwind Explicit";
        case FTCS:
            return "FTCS";
        case LAX_FRIEDRICHS:
            return "Lax-Friedrichs";
        case LAX_WENDROFF:
            return "Lax-Wendroff";
        case LEAPFROG:
            return "Leap-frog";
        case MACCORMACK:
            return "MacCormack";
        case ROE:
            return "Roe";
        case BEAM_WARMING:
            return "Beam-Warming Explicit";
        case IMPLICIT_UPWIND:
            return "Upwind Implicit";
        case CRANK_NICOLSON:
            return "Crank-Nicolson";
        case BEAM_WARMING_IMP:
            return "Beam-Warming Implicit";
        case RK3_WENO5:
            return "RK3+WENO5";
        default:
            return "Unknown";
        }
    }

    // 初始化函数
    bool initialize()
    {
        // 计算网格步长
        dx = (x2 - x1) / (ni - 1);

        try
        {
            // 分配内存
            x = make_unique<double[]>(ni);
            u = make_unique<double[]>(ni);
            um = make_unique<double[]>(ni);
            un = make_unique<double[]>(ni);

            // 初始化网格
            for (int i = 0; i < ni; i++)
            {
                x[i] = x1 + i * dx;
            }

            return true;
        }
        catch (const bad_alloc &e)
        {
            cerr << "内存分配失败: " << e.what() << endl;
            return false;
        }
    }
}