// GraduateProject_DNASelfAssembling.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"

int programInitiation(){
	ornt2bpornt_init();
	unitcell_init();
	return 0;
}
int _tmain(int argc, _TCHAR* argv[])
{
	programInitiation();
	patchPreset();


	random_device rd;
	mt19937 gen(123456);
	discrete_distribution<> d({ 40, 10, 10, 40 });
	map<int, int> m;
	for (int n = 0; n<10000; ++n) {
		++m[d(gen)];
	}
	for (auto p : m) {
		cout << p.first << " generated " << p.second << " times\n";
	}


	return 0;
}

