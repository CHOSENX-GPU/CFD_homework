///////////////////////////////////////////////////////////////////
// Output.cpp - 结果输出与可视化
///////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>
#include "Module.h"

using namespace CFDParams;
using namespace std;
namespace fs = filesystem;

// 结果输出函数
void call_output(double num, int idx, double cfl) {
    // 获取当前使用的求解方法名称
    string methodName = getSolverName(currentMethod);
    
    // 构建输出目录结构
    string outerFolderName = "./output" + to_string(idx) + "_" + methodName;
    ostringstream oss;
    oss << fixed << setprecision(2) << cfl;
    string cfl_format = oss.str();
    string innerFolderName = outerFolderName + "/CFL_" + cfl_format;
    
    try {
        // 创建外层目录
        if (!fs::exists(outerFolderName)) {
            fs::create_directory(outerFolderName);
        }
        
        // 创建内层目录
        if (!fs::exists(innerFolderName)) {
            fs::create_directory(innerFolderName);
        }
        
        // 构建文件名
        string filename = innerFolderName + "/file" + to_string(static_cast<int>(num)) + ".dat";
        
        // 打开输出文件
        ofstream outfile(filename);
        if (!outfile) {
            cerr << "无法创建输出文件: " << filename << endl;
            return;
        }
        
        // 写入Tecplot格式头
        outfile << "VARIABLEs = X,u" << '\n';
        outfile << "ZONE I=" << ni << '\n';
        outfile << "datapacking=block" << '\n';
        
        // 写入坐标数据
        for (int i = 0; i < ni; i++) {
            outfile << setprecision(8) << x[i] << '\n';
        }
        
        // 写入速度数据
        for (int i = 0; i < ni; i++) {
            outfile << setprecision(8) << u[i] << '\n';
        }
        
        outfile.close();
        cout << "输出结果已保存至: " << filename << endl;
    }
    catch (const exception& e) {
        cerr << "输出文件时发生错误: " << e.what() << endl;
    }
}