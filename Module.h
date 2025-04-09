//////////////////////////////////////////////////////////////////////////////////////
// Module.h，用于声明模块的全局变量
/////////////////////////////////////////////////////////////////////////////

#pragma once

extern int ni; //网格节点数
extern double x1; //起始坐标
extern double x2; //终止坐标
extern double dx; //计算网格步长
extern double *x; //网格节点坐标

extern double *u; //存储n时刻的速度
extern double *um; //过渡项，存储n+1时刻的速度
extern double *un; //存储n-1时刻的速度

extern double t0; //起始时间
extern double tout; //终止时间
extern double dt; //计算时间步长

extern double cfl; //CFL数
extern double a; //波速

extern int nloop; //循环步数
extern int nout; //输出间隔
