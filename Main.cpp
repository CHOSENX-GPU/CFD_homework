#include<iostream>
#include"Module.h"
#include"Main.h"

using namespace std;

int main()
{
    int idx=10;//定义输出名

    call_mesh1d(); //调用网格生成函数
    call_init(); //调用初始化函数
    call_output(0,idx,cfl);// 输出初始化结果检查是否正确
    call_CFL(); //调用CFL条件函数
    do
    {
        call_boundary(); //调用边界条件函数
        call_solve_10();//调用求解函数
        t0+=dt; //更新时间
        nloop = nloop + 1;
        if (nloop % nout == 0) //每隔nout步输出一次结果
		{
			call_output(nloop,idx,cfl);
		}
    } while (t0<=tout);

    // call_release(); //释放空间

    return 0;
}