//////////////////////////////////////////////////////////////////////////
// Module.cpp，用于定义模块的全局变量
//////////////////////////////////////////////////////////////////////////

int ni=101; //网格节点数
double x1=0.0; //起始坐标
double x2=1.0; //终止坐标
double dx; //计算网格步长

double *x=nullptr; //网格节点坐标    
double *u=nullptr; //存储n时刻的速度
double *um=nullptr; //过渡项，存储n+1时刻的速度
double *un=nullptr; //存储n-1时刻的速度

double t0=0.0; //起始时间
double tout=3; //终止时间
double dt; //计算时间步长

double cfl=0.8; //CFL数
double a=0.2; //波速

int nloop = 0; //循环步数
int nout = 2; //输出间隔

