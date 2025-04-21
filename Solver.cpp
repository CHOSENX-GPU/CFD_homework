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

// double TDMA(double *a, double *b, double *c, double *d, int n) // Thomas算法
// {
//     double *x = new double[n]; // 存储解的数组
//     double *p = new double[n]; // 存储中间变量c1
//     double *q = new double[n]; // 存储中间变量d1

//     p[0] = c[0] / b[0];
//     q[0] = d[0] / b[0];

//     for (int i = 1; i < n; i++)
//     {
//         double denom = b[i] - a[i] * p[i - 1];
//         p[i] = c[i] / denom;
//         q[i] = (d[i] + a[i] * q[i - 1]) / denom;
//     }

//     x[n - 1] = q[n - 1];
//     for (int i = n - 2; i >= 0; i--)
//     {
//         x[i] = p[i] * x[i + 1] + q[i];
//     }

//     delete[] p;
//     delete[] q;

//     return *x;
// }

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
    for (int i = 1; i < ni - 1; i++)
    {
        um[i] = u[i] - a * dt / (2 * dx) * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni - 1; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_3() // Lax-Friedrichs格式
{
    for (int i = 1; i < ni - 1; i++)
    {
        um[i] = 0.5 * (u[i + 1] + u[i - 1]) - 0.5 * a * dt / dx * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni - 1; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_4() // Lax-Wendroff格式
{
    for (int i = 1; i < ni - 1; i++)
    {
        um[i] = u[i] - 0.5 * a * dt / dx * (u[i + 1] - u[i - 1]) + 0.5 * a * a * dt * dt / (dx * dx) * (u[i + 1] - 2 * u[i] + u[i - 1]); // 递推计算n+1时刻的速度
    }

    for (int i = 1; i < ni - 1; i++)
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
        for (int i = 1; i < ni - 1; i++)
        {
            um[i] = un[i] - a * dt / dx * (u[i + 1] - u[i - 1]); // 递推计算n+1时刻的速度
        }
    }
    for (int i = 1; i < ni - 1; i++)
    {
        un[i] = u[i]; // 更新速度
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_6() // MacCormack格式
{
    for (int i = 1; i < ni - 1; i++) // 校正步
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
    for (int i = 1; i < ni - 1; i++)
    {
        fr = 0.5 * a * (u[i] + u[i + 1]) - 0.5 * fabs(a) * (u[i + 1] - u[i]);
        fl = 0.5 * a * (u[i] + u[i - 1]) - 0.5 * fabs(a) * (u[i] - u[i - 1]);
        um[i] = u[i] - dt / dx * (fr - fl); // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni - 1; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_8() // Beam-Warming显示格式
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

void call_solve_10() // Crank-Nicolson格式
{
    double a1 = a * dt / (2 * dx);
    double a2 = 1.0 + a1;
    um[0] = u[0]; // 处理第一个点
    for (int i = 1; i < ni; i++)
    {
        um[i] = (u[i] + a1 * (um[i - 1] + u[i - 1] - u[i])) / a2; // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

void call_solve_11() // Beam-Warming隐式格式
{
    um[0] = u[0];                                                                           // 处理第一个点
    um[1] = (u[1] + a * dt / (2 * dx) * (um[0] + u[0] - u[1])) / (1.0 + a * dt / (2 * dx)); // 处理第二个点
    double a1 = -a * dt / (2 * dx);
    double a2 = 2 * a * dt / dx;
    double a3 = 1.0 + a1 + a2;
    for (int i = 2; i < ni; i++)
    {
        um[i] = (u[i] + a1 * um[i - 2] + a2 * um[i - 1]) / a3; // 递推计算n+1时刻的速度
    }
    for (int i = 1; i < ni; i++)
    {
        u[i] = um[i]; // 更新速度
        // cout<<i<<'\t'<<u[i]<<endl;
    }
    // system("pause");
}

// 五阶WENO格式
double WENO(double vmm, double vm, double v, double vp, double vpp)
{
    const double eps = 1e-6;
    double S1 = 13.0 / 12.0 * (vmm - 2.0 * vm + v) * (vmm - 2.0 * vm + v) + 3.0 / 12.0 * (vmm - 4.0 * vm + 3.0 * v) * (vmm - 4.0 * vm + 3.0 * v);
    double S2 = 13.0 / 12.0 * (vm - 2.0 * v + vp) * (vm - 2.0 * v + vp) + 3.0 / 12.0 * (vm - vp) * (vm - vp);
    double S3 = 13.0 / 12.0 * (v - 2.0 * vp + vpp) * (v - 2.0 * vp + vpp) + 3.0 / 12.0 * (3.0 * v - 4.0 * vp + vpp) * (3.0 * v - 4.0 * vp + vpp);

    double alpha0 = 0.1 / ((eps + S1) * (eps + S1));
    double alpha1 = 0.6 / ((eps + S2) * (eps + S2));
    double alpha2 = 0.3 / ((eps + S3) * (eps + S3));

    double w0 = alpha0 / (alpha0 + alpha1 + alpha2);
    double w1 = alpha1 / (alpha0 + alpha1 + alpha2);
    double w2 = alpha2 / (alpha0 + alpha1 + alpha2);

    double u0 = (1.0 / 3.0) * vmm - (7.0 / 6.0) * vm + (11.0 / 6.0) * v;
    double u1 = -(1.0 / 6.0) * vm + (5.0 / 6.0) * v + (1.0 / 3.0) * vp;
    double u2 = (1.0 / 3.0) * v + (5.0 / 6.0) * vp - (1.0 / 6.0) * vpp;

    return w0 * u0 + w1 * u1 + w2 * u2;
}

// Ghost点处理函数
double bound(int i, double *u)
{
    if (i < 0)
    {
        return u[0];
    }
    else if (i > ni - 1)
    {
        return u[ni - 1];
    }
    else
    {
        return u[i];
    }
}

// 计算数值通量
double flux(int i, double *u)
{
    double fmm, fm, f, fp, fpp;
    if (a < 0)
    {
        fmm = (bound(i + 3, u) - bound(i + 2, u)) / dx;
        fm = (bound(i + 2, u) - bound(i + 1, u)) / dx;
        f = (bound(i + 1, u) - bound(i, u)) / dx;
        fp = (bound(i, u) - bound(i - 1, u)) / dx;
        fpp = (bound(i - 1, u) - bound(i - 2, u)) / dx;
    }
    else
    {
        fmm = (bound(i - 2, u) - bound(i - 3, u)) / dx;
        fm = (bound(i - 1, u) - bound(i - 2, u)) / dx;
        f = (bound(i, u) - bound(i - 1, u)) / dx;
        fp = (bound(i + 1, u) - bound(i, u)) / dx;
        fpp = (bound(i + 2, u) - bound(i + 1, u)) / dx;
    }
    return -a * WENO(fmm, fm, f, fp, fpp);
}

// 3阶 Runge-Kutta法 + WENO5格式
void call_solve_12()
{
    double *u1 = new double[ni]; // 存储第一步结果
    double *u2 = new double[ni]; // 存储第二步结果
    // 第一步
    for (int i = 1; i < ni - 3; i++)
    {
        u1[i] = u[i] + dt * flux(i, u);
    }

    // 第二步
    for (int i = 1; i < ni - 1; i++)
    {

        u2[i] = 0.75 * u[i] + 0.25 * u1[i] + 0.25 * dt * flux(i, u1);
    }

    // 第三步
    for (int i = 1; i < ni - 1; i++)
    {

        um[i] = 1.0 / 3.0 * u[i] + 2.0 / 3.0 * u2[i] + 2.0 / 3.0 * dt * flux(i, u2);
    }

    // 更新速度
    for (int i = 1; i < ni - 1; i++)
    {

        u[i] = um[i];
    }
    // system("pause");
}