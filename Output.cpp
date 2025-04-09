///////////////////////////////////////////////////////////////////
//  输出结果
///////////////////////////////////////////////////////////////////

#include<iostream>
#include<iomanip>
#include<fstream>
#include<string>
#include<cstring>
#include<sstream> // 包含stringstream头文件
#include <filesystem>
#include"Module.h"

using namespace std;
namespace fs = filesystem;


void call_output(double num, int idx, double cfl)
{
    // 构建文件夹名称
    string outerFolderName = "./output" + to_string(idx);
    ostringstream oss;
    oss << fixed << setprecision(2) << cfl; // 设置小数点后两位精度
    string cfl_format = oss.str();
    string innerFolderName = outerFolderName + "/CFL_" + cfl_format;

    // 检查外层文件夹是否存在，不存在则创建
    if (!fs::exists(outerFolderName)) 
    {
        if (!fs::create_directory(outerFolderName)) 
        {
            cerr << "无法创建 " << outerFolderName << " 文件夹" << std::endl;
            return;
        }
    }

    // 检查内层文件夹是否存在，不存在则创建
    if (!fs::exists(innerFolderName)) 
    {
        if (!fs::create_directory(innerFolderName)) 
        {
            cerr << "无法创建 " << innerFolderName << " 文件夹" << std::endl;
            return;
        }
    }


    /*生成输出文件*/
    char filename[50] = { 0 }; char str[5];		/*filename--输出文件名，str--文件标识符*/
    string fullPath = innerFolderName + "/file";
    strcpy_s(filename, fullPath.c_str());
    _itoa_s(num, str, 10);						        /*num转字符串*/
    strcat_s(filename, str);			         		/*将str添加到filename结尾*/
    strcat_s(filename, ".dat");					      /*将''.dat''添加到filename结尾*/

    ofstream outfile(filename);
    outfile << "VARIABLEs = X,u" << '\n';
    outfile << "ZONE I=" << ni << '\n';
    outfile << "datapacking=block" << '\n';
    for (int i = 0; i < ni; i++)
        outfile << x[i] << '\n';
    for (int i = 0; i < ni; i++)
        outfile << u[i] << '\n';
    outfile.close();
} 