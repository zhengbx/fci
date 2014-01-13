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

#include "fciopt.h"
#include "math.h"
#include <iostream>
#include <cmath>
#include <omp.h>
#include <mkl.h>
#include <cstring>
#include <fstream>
#include "hamiltonian.h"

Davidson_Updator::Davidson_Updator(int max_dim, int _nA, int _nB, double** init_x, double** init_y, double** _D, double* _oei, double* _tei, int _nmo, int* _llistA, int*** _listA, int* _llistB, int*** _listB): D(_D), nA(_nA), nB(_nB), nmo(_nmo), dim(max_dim), it(0), oei(_oei), tei(_tei), llistA(_llistA), llistB(_llistB), listA(_listA), listB(_listB), osh(llistB != NULL) {
  G = new double[ind1(dim,0)];
  x.push_back(init_x);
  y.push_back(init_y);
  res = new double*[nA];
  for (int i = 0; i < nA; i++) {
    res[i] = new double[nB];
  }
}

double Davidson_Updator::inner(double** a, double** b) {
  double sum = 0;
  #pragma omp parallel default(none) shared(a, b) reduction(+:sum)
  {
  #pragma omp for schedule(guided, 2)
  for (int i = 0; i < nA; i++) {
    for (int j = 0; j < nB; j++) {
      sum += a[i][j] * b[i][j];
    }
  }
  }
  return sum;
}

int Davidson_Updator::add(double** &x, double* const vec, std::vector <double**> const &basis, int n) {
  #pragma omp parallel default(none) shared(x, basis, n)
  {
  for (int i = 0; i < n; i++) {
    #pragma omp for schedule(guided, 2)
    for (int j = 0; j < nA; j++) {
      for (int k = 0; k < nB; k++) {
        x[j][k] += basis.at(i)[j][k] * vec[i];
      }
    }
  }
  }
  return 0;
}

int Davidson_Updator::correct(double** &delta, double** diag, double** residual, double E) {
  #pragma omp parallel default(none) shared(delta, diag, residual, E)
  {
  #pragma omp for schedule(guided, 2)
  for (int i = 0; i < nA; i++) {
    for (int j = 0; j < nB; j++) {
      delta[i][j] = residual[i][j] / (E-diag[i][j] + 0.2);
    }
  }
  }
  return 0;
}

int Davidson_Updator::update(double** &x1, double** &y1, double& E, double& r) {
  if (it >= dim) {
    return -1;
  }
  // update G
  for (int i = 0; i <= it; i++) {
    G[ind1(it, i)] = inner(y.back(), x.at(i));
  }
  it++;
  // diagonalize G, get lowest eigenvalue/vector
  double* workG = new double[it*it];
  for (int i = 0; i < it; i++) {
    for (int j = 0; j <= i; j++) {
      workG[i*it+j] = G[ind1(i, j)];
    }
  }
  int m;
  double *evalue = new double[it], *evectors = new double[it*it];
  int *ifail = new int[it];
  LAPACKE_dsyevx(LAPACK_ROW_MAJOR, 'V', 'I', 'L', it, workG, it, 0, 0, 1, 1, 1e-10, &m, evalue, evectors, it, ifail);
  // eigenvalue
  E = evalue[0];
  delete []workG;
  delete []evalue;
  double *evector = new double[it];
  for (int i = 0; i < it; i++) {
    evector[i] = evectors[i*it];
  }
  delete []evectors;
  delete []ifail;
  // calculate eigenvector (just for output, not stored here)
  for (int i = 0; i < nA; i++) {
    memset(x1[i], 0, sizeof(double) * nB);
    memset(y1[i], 0, sizeof(double) * nB);
  }
  add(x1, evector, x, it);
  add(y1, evector, y, it);
  delete []evector;
  // make residual vector
  double v[2] = {1.0, -E};
  std::vector <double **> u;
  u.push_back(y1);
  u.push_back(x1);
  for (int i = 0; i < nA; i++) {
    memset(res[i], 0, sizeof(double) * nB);
  }
  add(res, v, u, 2);
  u.clear();
  r = std::sqrt(inner(res, res));
  // construct correction vector
  double** x0 = new double*[nA];
  double** y0 = new double*[nA];
  for (int i = 0; i < nA; i++) {
    x0[i] = new double[nB];
    y0[i] = new double[nB];
  }
  correct(x0, D, res, E);
  normalize(x0, nA, nB);
  // Schmidt orthogonalization
  double *ips = new double[it];
  for (int i = 0; i < it; i++) {
    ips[i] = -inner(x0, x.at(i));
  }
  add(x0, ips, x, it);
  // normalization
  normalize(x0, nA, nB);
  x.push_back(x0);
  if (osh) {
    Hamp(y0, x0, D, oei, tei, nmo, nA, llistA, listA, nB, llistB, listB);
  } else {
    Hamp(y0, x0, D, oei, tei, nmo, nA, llistA, listA);
  }
  y.push_back(y0);
  return 0;
}

Davidson_Updator::~Davidson_Updator() {
  delete[] G;
  for (int i = 1; i < x.size(); i++) {
    for (int j = 0; j < nA; j++) {
      delete[] x.at(i)[j];
      delete[] y.at(i)[j];
    }
    delete[] x.at(i);
    delete[] y.at(i);
  }
  x.clear();
  y.clear();
  for (int i = 0; i < nA; i++) {
    delete[] res[i];
  }
  delete[] res;
}
