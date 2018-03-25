// stdafx.h : 标准系统包含文件的包含文件，
// 或是经常使用但不常更改的
// 特定于项目的包含文件
//

#pragma once

#include <ctime>
#include <random>

using namespace std;

// The following will be initialized at main
extern mt19937 gen; ///< Main random generator

extern clock_t t_start;
extern clock_t t_end;
