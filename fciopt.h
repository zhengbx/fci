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

#ifndef OPT
#define OPT
#include <vector>

class Davidson_Updator {
  double *G; // subspace Hamiltonian as triangular matrix
  std::vector <double**> x, y; // vectors
  double **res, **D; // correction vector
  int nA, nB, nmo;
  int dim, it;
  double *oei, *tei;
  int *llistA, *llistB, ***listA, ***listB;
  bool osh;
  double inner(double** a, double** b);
  int add(double** &x, double* const vec, std::vector <double**> const &basis, int n);
  int correct(double** &delta, double** diag, double** residual, double E);
public:
  Davidson_Updator(int max_dim, int _nA, int _nB, double** init_x, double** init_y, double** _D, double* _oei, double* _tei, int _nmo, int* _llistA, int*** _listA, int* _llistB, int*** _listB);
  int update(double** &x1, double** &y1, double& E, double& r);  // error if reach the max dimension
  ~Davidson_Updator();
};
#endif
