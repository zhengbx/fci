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

// calculate factorial
int fact(int n);
// calculation combination number C_m^n
int comb(int n, int m);
// convert occupation bit string into address
int add(int nbit, int str);
// convert address into occupation bit string
int str(int nbit, int nocc, int add);
// insert a bit to a string
int insbit(int pstr, int site, int occ);
// count occupied orbitals between two indicated sites (themselves not counted)
int countocc(int str, int i, int j);
// calculate time difference
double dtime(struct timeval t_start, struct timeval t_end);
// calculate index(i,j)=i*(i+1)/2+j (i >= j)
int ind1(int i, int j);
// calculate index(i,j)=i*(i-1)/2+j (i>j)
int ind2(int i, int j);
// normalize a vector
int normalize(double** &vec, int m, int n);
