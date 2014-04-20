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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdio>
#include <cmath>
#include "math.h"
#include "hamiltonian.h"
#include "timer.h"
#include "fciopt.h"

const int max_iter = 20;
const int max_dim = 6;
const double tolE = 1e-7;
const double tolr = 1e-4;

int readinput(char* filename, int& nelecA, int& nelecB, int& nmo, int& nA, int& nB, double* &tei, double* &oei, double& E0, int& ntei, int& noei, bool& open) {
  std::ifstream input(filename);
  int nelec, mz;
  input >> nmo >> nelec >> mz;
  if ((nelec + mz) % 2 != 0) {
    std::cout << "Error: Wrong spin configuration.\n";
    return -1;
  }
  if (mz == 0) {
    open = false;
  } else {
    open = true;
  }
  nelecA = (nelec + mz) / 2;
  nelecB = (nelec - mz) / 2;
  printf("Full CI Calculation\n\n%-20s=  %d    %-25s=  %d    %-25s=  %d\n", "no. of orbitals", nmo, "no. of alpha electrons", nelecA, "no. of beta electrons", nelecB);
  nA = comb(nmo, nelecA);
  nB = comb(nmo, nelecB);
  printf("Number of alpha strings   =  %d   Number of beta strings   =  %d\n", nA, nB);
  printf("Memory required for FCI vector  =  %10.4f MB\n", (double) (nA * nB) / 128 /1024);
  noei = ind1(nmo, 0);
  oei = new double[noei];
  memset(oei, 0, sizeof(double) * noei);
  ntei = ind1(ind1(nmo, 0), 0);
  tei = new double[ntei];
  memset(tei, 0, sizeof(double) * ntei);
  E0 = 0;
  double temp;  
  while (input >> temp) {
    int i, j, k, l;
    input >> i >> j >> k >> l;
    if (k != 0 && l != 0) {
      tei[ind1(ind1(i - 1, j - 1), ind1(k - 1, l - 1))] = temp;
    } else if (i != 0 && j != 0) {
      oei[ind1(i - 1, j - 1)] = temp;
    } else {
      E0 = temp;
    }
  }
  input.close();
  return 0;
}

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cout << "Input File Not Found\n";
    return -1;
  }
  Timer t_total("the whole computation");
  t_total.start();
  int nelecA, nelecB, nmo, ntei, noei;
  int nA, nB;
  double *tei, *oei;
  double E0;
  bool osh;
  // read the input file
  // get everything we need for the calculation, number of orbitals, electron occupation
  // integrals, nuclear energy, etc.
  if (readinput(argv[1], nelecA, nelecB, nmo, nA, nB, tei, oei, E0, ntei, noei, osh) == -1) {
    std::cout << "Input File Error\n";
    return -1;
  }
  // construct mapping lists
  int ***listA, ***listB = NULL; // lists O, OV, OO, OOV, OOVV
  int *llistA, *llistB = NULL;  // corresponding list lengths
  listA = new int**[5];
  llistA = new int[5];
  Timer t_list("making Hamiltonian mapping lists");
  t_list.start();
  if (!osh) {
    osh = true; // closed shell not implemented correctly
    nB = nA;
    nelecB = nelecA;
  }
  if (osh) {
    printf("\ndo open shell calculation...\n\n");
  }
  for (list_type i = 0; i < 5; i++) {
    makelist(listA[i], llistA[i], nmo, nelecA, i);
  }
  if (osh) {
    listB = new int**[5];
    llistB = new int[5];
    for (list_type i = 0; i < 5; i++) {
      makelist(listB[i], llistB[i], nmo, nelecB, i);
    }
  } else {
    printf("\ndo closed shell calculation...\n\n");
  }
  t_list.finish();
  t_list.print();  
  // make diagonal elements
  Timer t_diag("making diagonal Hamiltonian");
  t_diag.start();
  double** D;
  if (osh) {
    makeDiagH(D, nmo, oei, tei, nA, llistA, listA, nB, llistB, listB);
  } else {
    makeDiagH(D, nmo, oei, tei, nA, llistA, listA);
  }
  
  t_diag.finish();
  t_diag.print();
  // make initial guess vector x (x0)
  // perturbed HF vector
  Timer t_init("initializing guess vectors");
  t_init.start();
  double **x0, **y0 = NULL;
  if (osh) {
    init_X(x0, nA, nB);
    Hamp(y0, x0, D, oei, tei, nmo, nA, llistA, listA, nB, llistB, listB);
  } else {
    init_X(x0, nA);
    Hamp(y0, x0, D, oei, tei, nmo, nA, llistA, listA);
  }
  t_init.finish();
  t_init.print();
  // Davidson algorithm for ground state energy
  int iter = 0;
  if (!osh) {
    nB = nA;
    llistB = NULL;
    listB = NULL;
  }
  Davidson_Updator davidson(max_iter, nA, nB, x0, y0, D, oei, tei, nmo, llistA, listA, llistB, listB);
  // allocate memory for x, y
  double **x = new double*[nA], **y = new double*[nA], E_old = 0, E = 0;\
  for (int i = 0; i < nA; i++) {
    x[i] = new double[nB];
    y[i] = new double[nB];
  }
  bool converge = false;
  printf("\nIter.      Energy       Energe Change      Residual        Time\n");
  Timer t_iter("");

  for (iter = 0; iter < max_iter; iter++) {
    t_iter.start();
    double r = 1; // residual vector
    if (davidson.update(x, y, E, r) < 0) {
      std::cout << "Max iterative dimension reached.\n";
      break;
    }
    t_iter.finish();
    E += E0;
    printf("%3d%17.8f%15.2e%16.2e%12.2f\n", iter+1, E, E-E_old, r, t_iter.time());
    if (std::abs(E-E_old) < tolE && r < tolr) {
      converge = true;
      std::cout << "\nFCI Converged.\n";
      printf("E = %18.12f Hartree\n\n", E);
      break;
    }
    E_old = E;
  }
  
  if (true) {
    print_vec(x, nmo, nelecA, nelecB);
  }
  
  t_total.finish();
  t_total.print();
  // disposal the arrays
  delete[] oei;
  delete[] tei;
  return 0;
}
