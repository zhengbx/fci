/*
# This file is part of the Fci program. Fci is free software: you
# can redistribute it and/or modify it under the terms of the GNU
# General Public License as published by the Free Software Foundation,
# version 3.
# 
# Fci is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
# 
# Authors:
#    Bo-Xiao Zheng, 2013
*/

// calculate y=Hx

typedef int list_type;

const list_type
  O = 0,
  OV = 1,
  OO = 2,
  OOV = 3,
  OOVV = 4;

int makelist(int** &list, int& len, int nmo, int nocc, list_type type);

int makelistO(int* &plist, int i, int len, int nmo, int nocc);

int makelistOV(int* &plist, int i, int j, int len, int nmo, int nocc);

int makelistOO(int* &plist, int i, int j, int len, int nmo, int nocc);

int makelistOOV(int* &plist, int i, int j, int k, int len, int nmo, int nocc);

int makelistOOVV(int* &plist, int i, int j, int k, int l, int len, int nmo, int nocc);

int makeDiagH(double** &D, int nmo, double* oei, double* tei, int nA, int* llistA, int*** listA, int nB = 0, int* llistB = NULL, int*** listB = NULL);

int init_X(double** &x, int nA, int nB = 0);

int Hamp(double** &y, double** x, double** D, double* oei, double* tei, int nmo, int nA, int* llistA, int*** listA, int nB = 0, int* llistB = NULL, int*** listB = NULL);
