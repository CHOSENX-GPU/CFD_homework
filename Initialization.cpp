///////////////////////////////////////////////////////////////////
//  定义初始条件
///////////////////////////////////////////////////////////////////

#include<iostream>
#include<cmath>
#include"Module.h"

using namespace std;

// 定义初始条件(脉冲方波)
void call_init()
{
    u = new double[ni]; //分配内存，定义速度
    um = new double[ni]; //分配内存，定义n+1时刻的速度
    un = new double[ni]; //分配内存，定义n-1时刻的速度
    
    t0=0.0; //起始时间

    for (int i=0;i<ni;i++)
    {
        if (x[i]>=0.0 && x[i]<=0.3)
        {
            u[i]=1.0; //定义坐标小于0.3处的速度为1.0
        }
        else
        {
            u[i]=0.0; //其余坐标点处速度为0
        }
        // cout<<i<<'\t'<<u[i]<<endl;
    }
}