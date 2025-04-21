///////////////////////////////////////////////////////////////////
//  采用迎风格式求解一维线性对流方程
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "Module.h"

using namespace std;

double sign(double x)
{
    if (x > 0)
    {
        return 1.0;
    }
    else
    {
        return -1.0;
    }
}

void call_solve_1() // 时间向前空间向前差分格式（UpWind）
{
    for (int i = 1; i < ni; i++)
    {
        um[i] = u[i] - a * dt / dx * (u[i] - u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_2() // 时间向前空间中心差分格式（FTCS）
{
    for (int i = 1; i < ni; i++)
    {
        um[i] = u[i] - a * dt / (2 * dx) * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_3() // Lax-Friedrichs格式
{
    for (int i = 1; i < ni; i++)
    {
        um[i] = 0.5 * (u[i + 1] + u[i - 1]) - 0.5 * a * dt / dx * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_4() // Lax-Wendroff格式
{
    for (int i = 1; i < ni; i++)
    {
        um[i] = u[i] - 0.5 * a * dt / dx * (u[i + 1] - u[i - 1]) + 0.5 * a * a * dt * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_5() // Leap-frog格式
{
    if (nloop == 0) // 当nloop=0时，时间推进采用迎风格式
    {
        for (int i = 1; i < ni; i++)
        {
            um[i] = u[i] - a * dt / dx * (u[i] - u[i - 1]); // 递推计算n+1时刻的速度
        }
    }
    else // 当nloop>0时，um和un都为上一个时刻的值
    {
        for (int i = 1; i < ni; i++)
        {
            um[i] = un[i] - a * dt / dx * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
        }
    }
    for (int i = 1; i < ni; i++)
    {
        un[i] = u[i]; // 更新速度
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_6() // MacCormack格式
{
    for (int i = 1; i < ni; i++) // 校正步
    {
        um[i] = u[i] - a * dt / dx * (u[i + 1] - u[i]);
    }
    for (int i = 1; i < ni; i++) // 预测步
    {
        u[i] = 0.5 * (u[i] + um[i]) - 0.5 * a * dt / dx * (um[i] - um[i - 1]);
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_7() // Roe格式
{
    double fr, fl;
    for (int i = 1; i < ni; i++)
    {
        fr = 0.5 * a * (u[i] + u[i + 1]) - 0.5 * fabs(a) * (u[i + 1] - u[i]);
        fl = 0.5 * a * (u[i] + u[i - 1]) - 0.5 * fabs(a) * (u[i] - u[i - 1]);
        um[i] = u[i] - dt / dx * (fr - fl); // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_8() // Beam-Warming格式
{
    um[1] = u[1] - a * dt / dx * (u[1] - u[0]); // 处理第二个点
    for (int i = 2; i < ni; i++)
    {
        um[i] = u[i] - 0.5 * a * dt / dx * (3 * u[i] - 4 * u[i - 1] + u[i - 2]); // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_9() // 隐式迎风格式
{
    double a1 = a * dt / dx;
    double a2 = 1.0 + a1;
    um[0] = u[0]; // 处理第一个点
    for (int i = 1; i < ni; i++)
    {
        um[i] = (u[i] + a1 * um[i - 1]) / a2; // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_10() //Crank-Nicolson格式
{
    double a1 = a * dt / (2 * dx);
    double a2 = 1.0 + a1;
    um[0] = u[0]; // 处理第一个点
    for (int i = 1; i < ni; i++)
    {
        um[i] = (u[i] + a1 * (um[i - 1] + u[i - 1]- u[i])) / a2; // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}
