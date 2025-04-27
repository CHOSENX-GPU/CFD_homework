///////////////////////////////////////////////////////////////////
//  定义计算网格
///////////////////////////////////////////////////////////////////

#include<iostream>
#include"Module.h"

using namespace std;
using namespace CFDParams;

void call_mesh1d() {
    // 网格已在CFDParams::initialize()中生成，
    // 此函数保留以维持向后兼容性
    cout << "网格生成完成，dx = " << dx << endl;
}