// stdafx.h : ��׼ϵͳ�����ļ��İ����ļ���
// ���Ǿ���ʹ�õ��������ĵ�
// �ض�����Ŀ�İ����ļ�
//

#pragma once

#include <ctime>
#include <random>

using namespace std;

// The following will be initialized at main
extern mt19937 gen; ///< Main random generator

extern clock_t t_start;
extern clock_t t_end;
