void call_mesh1d(); // 网格划分函数
void call_init(); // 初始条件函数
void call_boundary(); // 边界条件函数
void call_CFL(); // CFL数计算函数
void call_solve_1(); // 一阶迎风格式
void call_solve_2(); // 时间向前空间中心差分格式（FTCS）
void call_solve_3(); // Lax-Friedrichs格式
void call_solve_4(); // Lax-Wendroff格式
void call_solve_5(); // Leap-frog格式
void call_solve_6(); // MacCormack格式
void call_solve_7(); // Beam-Warming显式格式
void call_solve_8(); // 隐式迎风格式
void call_solve_9(); // Crank-Nicolson格式
void call_solve_10(); // Beam-Warming隐式格式
void call_solve_11(); // 三阶Runge-Kutta + 五阶WENO格式
void callSolver(CFDParams::SolverMethod method); // 求解器选择函数
void call_output(double num, int idx, double cfl); // 结果输出函数