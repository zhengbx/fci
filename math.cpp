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

#include "math.h"
#include <iostream>
#include <cmath>
#include <omp.h>

using namespace std;

int fact(int n) {
  // n >= 0
  int f = 1;
  for (int i = n; i > 1; i--) {
    f *= i;
  }
  return f;
}

int comb(int n, int m) {
  // n >= m >= 0
  if (m > n || n < 0 || m < 0) {
    return 0;
  }
  if (m > n / 2) {
    return comb(n, n - m);
  }
  int c = 1;
  for (int i = n; i > n - m; i--) {
    c *= i;
  }
  c /= fact(m);
  return c;
}

int add(int nbit, int str) {
  int bit = 1;
  int addr = 0;
  int occ = 0;
  for (int k = 0; k < nbit; k++) {
    if ((bit & str) != 0) {
      occ += 1;
      addr += comb(k, occ);
    }
    bit *= 2;
  }
  return addr + 1;
}

int str(int nbit, int nocc, int add) {
  add--;
  int i = nocc;
  int string = 0;
  for (int j = nbit - 1; j >= 0; j--) {
    if (add >= comb(j, i)) {
      string += 1 << j;
      add -= comb(j, i);
      i--;
    }
  }
  return string;
}

int insbit(int pstr, int site, int occ) {
  // pstr: original partial string
  // site: site to insert
  // occupation at the site
  int sitelabel = 1 << site;
  return (pstr / sitelabel) * sitelabel * 2 + sitelabel * occ + pstr % sitelabel;
}

int countocc(int str, int i, int j) {
  if (i < j) {
    int temp;
    temp = i;
    i = j;
    j = temp;
  }
  // first extract the interesting part
  str = (str % (1 << i)) / (1 << (j + 1));
  // then count site by site
  int count = 0;
  int sitelabel = 1 << (i - j - 1);
  while (sitelabel > 0) {
    count += str / sitelabel;
    str -= (str / sitelabel) * sitelabel;
    sitelabel /= 2;
  }
  return count;
}

double dtime(struct timeval t_start, struct timeval t_end) {
  return (double)(t_end.tv_sec-t_start.tv_sec)+((double)(t_end.tv_usec-t_start.tv_usec))/1000000.;
}

int ind1(int i, int j) {
  return i>j ? i*(i+1)/2 + j : j*(j+1)/2 + i;
}

int ind2(int i, int j) {
  return i>j ? i*(i-1)/2 + j : j*(j-1)/2 + i;
}

int normalize(double** &vec, int m, int n) {
  double sum = 0;
  #pragma omp parallel default(none) shared(vec, m, n) reduction(+: sum)
  {
  #pragma omp for schedule(guided, 1)
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
	  sum += vec[i][j] * vec[i][j];
	}
  }
  }
  sum = sqrt(sum);
  #pragma omp parallel default(none) shared(vec, m, n, sum)
  {
  #pragma omp for schedule(guided, 1)
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      vec[i][j] /= sum;
	}
  }
  }
  return sum;
}
    
