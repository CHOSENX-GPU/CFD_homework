///////////////////////////////////////////////////////////////////
// Initialization.cpp - 初始条件设置
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include "Module.h"

using namespace CFDParams;
using namespace std;

// 设置初始条件(脉冲方波)
void call_init() {
    t0 = 0.0; // 初始时间
    
    try {
        // 初始速度分布(方波脉冲)
        for (int i = 0; i < ni; i++) {
            // 定义脉冲波：在[0.0, 0.3]区间内速度为1.0，其他位置为0
            u[i] = (x[i] >= 0.0 && x[i] <= 0.3) ? 1.0 : 0.0;
            // 初始化其他时间层
            um[i] = u[i];
            un[i] = u[i];
        }
        
        cout << "初始条件设置完成" << endl;
    }
    catch (const std::exception& e) {
        cerr << "初始化失败: " << e.what() << endl;
        exit(1);
    }
}