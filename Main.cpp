///////////////////////////////////////////////////////////////////
// Main.cpp - 程序主入口
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <chrono>
#include <iomanip>
#include "Module.h"
#include "Main.h"

using namespace CFDParams;
using namespace std;

int main() {
    cout << "===========================================" << endl;
    cout << "一维线性对流方程求解器" << endl;
    cout << "===========================================" << endl;
    
    // 设置求解方法 (1-12)
    int idx = 1;  // 迎风格式
    
    // 根据idx设置当前求解方法
    currentMethod = static_cast<SolverMethod>(idx);
    
    // 初始化参数
    if (!initialize()) {
        cerr << "初始化失败，程序终止" << endl;
        return 1;
    }
    
    // 输出当前使用的求解方法
    cout << "使用求解方法: " << getSolverName(currentMethod) << endl;
    
    // 计时开始
    auto start_time = chrono::high_resolution_clock::now();
    
    // 创建网格
    call_mesh1d();
    
    // 设置初始条件
    call_init();
    
    // 输出初始条件
    call_output(0, idx, cfl);
    
    // 计算时间步长
    call_CFL();
    
    cout << "开始迭代计算..." << endl;
    
    // 主计算循环
    do {
        // 应用边界条件
        call_boundary();
        
        // 数值求解
        callSolver(currentMethod);
        
        // 更新时间
        t0 += dt;
        nloop++;
        
        // 输出结果
        if (nloop % nout == 0) {
            cout << "迭代次数: " << nloop << ", 当前时间: " << t0 << endl;
            call_output(nloop, idx, cfl);
        }
        
    } while (t0 <= tout);
    
    // 计时结束
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    
    cout << "===========================================" << endl;
    cout << "计算完成！" << endl;
    cout << "总迭代次数: " << nloop << endl;
    cout << "模拟总时间: " << t0 << endl;
    cout << "计算耗时: " << duration.count() / 1000.0 << " 秒" << endl;
    cout << "===========================================" << endl;
    
    return 0;
}