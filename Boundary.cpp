///////////////////////////////////////////////////////////////////
//  定义边界条件
///////////////////////////////////////////////////////////////////

#include<iostream>
#include"Module.h"

using namespace std;

void call_boundary()
{
    u[0]=0.0; //左边界处速度为0
    u[ni]=0.0; //右边界处速度为0
}