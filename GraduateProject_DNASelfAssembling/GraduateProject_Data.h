#pragma once

extern const double NN_dH_kcal[4][4][4];
extern const double NN_dS_cal[4][4][4];
/*
	NN_dH_kcal[a][b][c], NN_dS_cal[a][b][c]

	5'  ---  a  ---  b  ---  3'
	         |       |
	3'  ---  a' ---  c  ---  5'

	where a' is the paired nucleotide of a.
	If c = b' then the index [a][b][b'] and [b'][a'][a] must refer to the same value.
*/
extern const double NN_dH_init_kcal;
extern const double NN_dS_init_cal;