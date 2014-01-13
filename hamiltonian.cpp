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

#include <cstdlib>
#include "hamiltonian.h"
#include "math.h"
#include <omp.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <assert.h>
#include <fstream>
#include "timer.h"

int makelist(int** &list, int& len, int nmo, int nocc, list_type type) {
  if (type == O) {
    list = new int*[nmo];
    len = comb(nmo - 1, nocc - 1);
    #pragma omp parallel default(none) shared(list, len, nmo, nocc)
    {
    #pragma omp for schedule(guided, 2)
    for (int i = 0; i < nmo; i++) {
      makelistO(list[i], i, len, nmo, nocc);
    }
    }
  } else if (type == OV) {
    list = new int*[nmo * nmo];
    len = comb(nmo - 2, nocc - 1);
    #pragma omp parallel default(none) shared(list, len, nmo, nocc)
    {
    #pragma omp for schedule(guided, 2)
    for (int m = 0; m < nmo * nmo; m++) {
      int i = m / nmo;
      int j = m % nmo;
      if (i != j) {
        makelistOV(list[m], i, j, len, nmo, nocc);
      }
    }
    }
  } else if (type == OO) {
    list = new int*[ind2(nmo, 0)];
    len = comb(nmo - 2, nocc - 2);
    #pragma omp parallel default(none) shared(list, len, nmo, nocc)
    {
    #pragma omp for schedule(guided, 1)
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        makelistOO(list[ind2(i,j)], i, j, len, nmo, nocc);
      }
    }
    }
  } else if (type == OOV) {
    list = new int*[nmo*nmo*nmo]; // i, j, k go over all different orbitals
    len = comb(nmo - 3, nocc - 2);
    #pragma omp parallel default(none) shared(list, len, nmo, nocc)
    {
    #pragma omp for schedule(guided, 1)
    for (int m = 0; m < nmo*nmo*nmo; m++) {
      int i = m / (nmo * nmo);
      int j = (m % (nmo * nmo)) / nmo;
      int k = m % nmo;
      if (i != j && i != k && j != k) {
        makelistOOV(list[m], i, j, k, len, nmo, nocc);
      }
    }
    }
  } else if (type == OOVV) {
    list = new int*[ind2(nmo, 0) * ind2(nmo, 0)];
    len = comb(nmo - 4, nocc - 2);
    #pragma omp parallel default(none) shared(list, len, nmo, nocc)
    {
    #pragma omp for schedule(guided, 1)
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        for (int k = 0; k < nmo; k++) {
          if (k == i || k == j) {
            continue;
          }
          for (int l = 0; l < k; l++) {
            if (l == i || l == j) {
              continue;
            }
            makelistOOVV(list[ind2(i,j)*ind2(nmo, 0)+ind2(k,l)], i, j, k, l, len, nmo, nocc);
          }
        }
      }
    }
    }
  }
  return 0;
}

int makelistO(int* &plist, int i, int len, int nmo, int nocc) {
  plist = new int[len];
  for (int s = 0; s < len; s++) {
    int pstr = str(nmo - 1, nocc - 1, s + 1);
    plist[s] = add(nmo, insbit(pstr, i, 1));
  }
  return 0;
}

int makelistOV(int* &plist, int i, int j, int len, int nmo, int nocc) {
  plist = new int[len];
  for (int s = 0; s < len; s++) {
    int pstr = str(nmo - 2, nocc - 1, s + 1);
    if (j < i) {
      pstr = insbit(pstr, j, 0);
      pstr = insbit(pstr, i, 1);
    } else {
      pstr = insbit(pstr, i, 1);
      pstr = insbit(pstr, j, 0);
    }
    int sign = 1 - (countocc(pstr, i, j) % 2) * 2;
    plist[s] = sign * add(nmo ,pstr);
  }
  return 0;
}

int makelistOO(int* &plist, int i, int j, int len, int nmo, int nocc) {
  plist = new int[len];
  for (int s = 0; s < len; s++) {
    int pstr = str(nmo - 2, nocc - 2, s + 1);
    pstr = insbit(pstr, j, 1);
    pstr = insbit(pstr, i, 1);
    plist[s] = add(nmo, pstr);
    //cout << OO[ind2(i, j)][s] << '\t';
  }
  return 0;
}

int makelistOOV(int* &plist, int i, int j, int k, int len, int nmo, int nocc) {
  plist = new int[len];
  for (int s = 0; s < len; s++) {
    int pstr = str(nmo - 3, nocc - 2, s + 1);
    if (i < j && j < k) {
      pstr = insbit(pstr, i, 1);
      pstr = insbit(pstr, j, 1);
      pstr = insbit(pstr, k, 0);
    } else if (i < k && k < j) {
      pstr = insbit(pstr, i, 1);
      pstr = insbit(pstr, k, 0);
      pstr = insbit(pstr, j, 1);
    } else if (j < i && i < k) {
      pstr = insbit(pstr, j, 1);
      pstr = insbit(pstr, i, 1);
      pstr = insbit(pstr, k, 0);
    } else if (j < k && k < i) {
      pstr = insbit(pstr, j, 1);
      pstr = insbit(pstr, k, 0);
      pstr = insbit(pstr, i, 1);
    } else if (k < i && i < j) {
      pstr = insbit(pstr, k, 0);
      pstr = insbit(pstr, i, 1);
      pstr = insbit(pstr, j, 1);
    } else {
      pstr = insbit(pstr, k, 0);
      pstr = insbit(pstr, j, 1);
      pstr = insbit(pstr, i, 1);
    }
    int sign = 1 - (countocc(pstr, j, k) % 2) * 2;
    plist[s] = add(nmo, pstr) * sign;
  }
  return 0;
}

int makelistOOVV(int* &plist, int i, int j, int k, int l, int len, int nmo, int nocc) {
  plist = new int[len];
  for (int s = 0; s < len; s++) {
    int pstr = str(nmo - 4, nocc - 2, s + 1);
    int sign = 0;
    if (l < j) {
      pstr = insbit(pstr, l, 0);
      if (k < j) {
        pstr = insbit(pstr, k, 0);
        pstr = insbit(pstr, j, 1);
        pstr = insbit(pstr, i, 1);
        sign = 1; // k between l, j
      } else {
        pstr = insbit(pstr, j, 1);
        if (k < i) {
          pstr = insbit(pstr, k, 0);
          pstr = insbit(pstr, i, 1);
        } else {
          pstr = insbit(pstr, i, 1);
          pstr = insbit(pstr, k, 0);
        }
      }
    } else {
      pstr = insbit(pstr, j, 1);
      if (i < l) {
        pstr = insbit(pstr, i, 1);
        pstr = insbit(pstr, l, 0);
        pstr = insbit(pstr, k, 0);
        sign = 1; // i between l, j
      } else {
        pstr = insbit(pstr, l, 0);
        if (k < i) {
          pstr = insbit(pstr, k, 0);
          pstr = insbit(pstr, i, 1);
        } else {
          pstr = insbit(pstr, i, 1);
          pstr = insbit(pstr, k, 0);
        }
      }
    }
    sign += countocc(pstr, i, k) + countocc(pstr, j, l);
    sign = 1 - (sign % 2) * 2;
    plist[s] = add(nmo, pstr) * sign;
  }
  return 0;
}

int makeDiagH(double** &D, int nmo, double* oei, double* tei, int nA, int* llistA, int*** listA, int nB, int* llistB, int*** listB) {
  D = new double*[nA];
  if (nB == 0) { // closed shell
    // only half of the coefficients will be used
    for (int i = 0; i < nA; i++) {
      D[i] = new double[nA];
      memset(D[i], 0, sizeof(double) * nB);
    }
  } else {  // open shell
    // allocate memory
    for (int i = 0; i < nA; i++) {
      D[i] = new double[nB];
      memset(D[i], 0, sizeof(double) * nB);
    }
    #pragma omp parallel default(none) shared(nmo, oei, tei, llistA, llistB, listA, listB, nA, nB, D)
    {
    // use O(i) lists    
    // 1e operator, case i = j
    for (int i = 0; i < nmo; i++) {
      double hii = oei[ind1(i, i)];
      #pragma omp for schedule(guided, 2)
      for (int sA = 0; sA < llistA[O]; sA++) {
        int addrA = listA[O][i][sA] - 1;
        for (int addrB = 0; addrB < nB; addrB++) {
          D[addrA][addrB] += hii;
        }
      }
      #pragma omp barrier
      #pragma omp for schedule(guided, 2)
      for (int sB = 0; sB < llistB[O]; sB++) {
        int addrB = listB[O][i][sB] - 1;
        for (int addrA = 0; addrA < nA; addrA++) {
          D[addrA][addrB] += hii;
        }
      }
    }
    // 2e operator, alpha-beta part
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < nmo; j++) {
        double Jij = tei[ind1(ind1(i, i), ind1(j, j))];
        #pragma omp for schedule(guided, 2)
        for (int sA = 0; sA < llistA[O]; sA++) {
          int addrA = listA[O][i][sA] - 1;
          for (int sB = 0; sB < llistB[O]; sB++) {
            int addrB = listB[O][j][sB] - 1;
            D[addrA][addrB] += Jij;
          }
        }
      }
    }
    }
    // use OO(i, j) lists, 2e alpha-alpha beta-beta part
    #pragma omp parallel default(none) shared(nmo, tei, llistA, listA, llistB, listB, nA, nB, D) 
    {
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        double V = tei[ind1(ind1(i,i),ind1(j,j))] - tei[ind1(ind1(i, j), ind1(j, i))];
        #pragma omp for schedule(guided, 1)
        for (int sA = 0; sA < llistA[OO]; sA++) {
          int addrA = listA[OO][ind2(i,j)][sA] - 1;
          for (int addrB = 0; addrB < nB; addrB++) {
            D[addrA][addrB] += V;
          }
        }
        #pragma omp for schedule(guided, 1)
        for (int sB = 0; sB < llistB[OO]; sB++) {
          int addrB = listB[OO][ind2(i,j)][sB] - 1;
          for (int addrA = 0; addrA < nA; addrA++) {
            D[addrA][addrB] += V;
          }
        }
      }
    }
    }
  }
  return 0; 
}

int init_X(double** &x, int nA, int nB) {
  x = new double*[nA];
  if (nB == 0) {
    nB = nA;
  }
  for (int i = 0; i < nA; i++) {
    x[i] = new double[nB];
  }
  srand(time(NULL));
  for (int i = 0; i < nA; i++) {
    for (int j = 0; j < nB; j++) {
      x[i][j] = 0.2 * (((double) rand()) / RAND_MAX - 0.5) / nA / nB;
    }
  }
  x[0][0] = 1;
  normalize(x, nA, nB);
  return 0;
}

int Hamp(double** &y, double** x, double** D, double* oei, double* tei, int nmo, int nA, int* llistA, int*** listA, int nB, int* llistB, int*** listB) {
  bool osh = (nB > 0);
  if (!osh) {
    nB = nA;
  }
  // initialize y
  if (y == NULL) {
    y = new double*[nA];
    for (int i = 0; i < nA; i++) {
      y[i] = new double[nB];
    }
  }
  if (osh) {  // open shell
    #pragma omp parallel default(none) shared(x, y, D, nA, nB, nmo, llistA, listA, llistB, listB, oei, tei)
    {
    // apply diagonal part
    Timer t_diag("diag");
    t_diag.start();
    #pragma omp for schedule(guided, 1)
    for (int i = 0; i < nA; i++) {
      for (int j = 0; j < nB; j++) {
        y[i][j] = x[i][j] * D[i][j];
      }
    }
    t_diag.finish();
    //t_diag.print();
    // now off-diagonal part
    // one-electron excitation
    Timer t_1eA("1eA");
    t_1eA.start();
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        double hij = oei[ind1(i, j)];
        #pragma omp for schedule(guided, 1)        
        for (int sA = 0; sA < llistA[OV]; sA++) {
          int addrA1 = abs(listA[OV][i*nmo+j][sA]) - 1;
          int sign = 2 * (listA[OV][i*nmo+j][sA] > 0) - 1;
          int addrA2 = abs(listA[OV][j*nmo+i][sA]) - 1;
          double hij_signed = hij * sign;
          for (int addrB = 0; addrB < nB; addrB++) {
            y[addrA1][addrB] += x[addrA2][addrB] * hij_signed;
            y[addrA2][addrB] += x[addrA1][addrB] * hij_signed;
          }
        }
      }
    }
    t_1eA.finish();
    //t_1eA.print();
    Timer t_1eB("1eB");
    t_1eB.start();
    #pragma omp for schedule(guided, 2)
    for (int addrA = 0; addrA < nA; addrA++) {
      for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < i; j++) {
          double hij = oei[ind1(i, j)];
          for (int sB = 0; sB < llistB[OV]; sB++) {
            int addrB1 = abs(listB[OV][i*nmo+j][sB]) - 1;
            int sign = 2 * (listB[OV][i*nmo+j][sB] > 0) - 1;
            int addrB2 = abs(listB[OV][j*nmo+i][sB]) - 1;
            double hij_signed = hij * sign;
            y[addrA][addrB1] += x[addrA][addrB2] * hij_signed;
            y[addrA][addrB2] += x[addrA][addrB1] * hij_signed;
          }
        }
      }
    }
    t_1eB.finish();
    //t_1eB.print();
    // 2e A^AB part
    Timer t_AA("2e A^AB_A");
    t_AA.start();
    for (int i = 0; i < nmo; i++) {
      #pragma omp for schedule(guided, 1)      
      for (int sA = 0; sA < llistA[O]; sA++) {
        int addrA = listA[O][i][sA] - 1;
        for (int j = 0; j < nmo; j++) {
          for (int k = 0; k < j; k++) {
            double V = tei[ind1(ind1(j, k), ind1(i, i))]; // (jk||ii)
            for (int sB = 0; sB < llistB[OV]; sB++) {
              int addrB1 = abs(listB[OV][j*nmo+k][sB]) - 1;
              int sign = 2 * (listB[OV][j*nmo+k][sB] > 0) - 1;
              int addrB2 = abs(listB[OV][k*nmo+j][sB]) - 1;
              y[addrA][addrB1] += x[addrA][addrB2] * V * sign;
              y[addrA][addrB2] += x[addrA][addrB1] * V * sign;
            }
          }
        }
      }
    }
    t_AA.finish();
    //t_AA.print();
    Timer t_AB("2e A^AB_B");
    t_AB.start();
    for (int j = 0; j < nmo; j++) {
      for (int k = 0; k < j; k++) {
        #pragma omp for schedule(guided, 1)
        for (int sA = 0; sA < llistA[OV]; sA++) {
          int addrA1 = abs(listA[OV][j*nmo+k][sA]) - 1;
          int sign = 2 * (listA[OV][j*nmo+k][sA] > 0) - 1;
          int addrA2 = abs(listA[OV][k*nmo+j][sA]) - 1;
          for (int i = 0; i < nmo; i++) {
            double V = tei[ind1(ind1(j, k), ind1(i, i))]; // (jk||ii)
            for (int sB = 0; sB < llistB[O]; sB++) {
              int addrB = listB[O][i][sB] - 1;
              y[addrA1][addrB] += x[addrA2][addrB] * V * sign;
              y[addrA2][addrB] += x[addrA1][addrB] * V * sign;
            }
          }
        }
      }
    }
    t_AB.finish();
    //t_AB.print();
    // 2e B^AB part
    Timer t_B("2e B^AB");
    t_B.start();
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        #pragma omp for schedule(guided, 1)
        for (int sA = 0; sA < llistA[OV]; sA++) {
          int addrA1 = abs(listA[OV][i*nmo+j][sA]) - 1;
          int signA = 2 * (listA[OV][i*nmo+j][sA] > 0) - 1;
          int addrA2 = abs(listA[OV][j*nmo+i][sA]) - 1;
          int ij = ind1(i, j);
          for (int k = 0; k < nmo; k++) {
            for (int l = 0; l < k; l++) {
              double V_A = tei[ind1(ij, ind1(k, l))] * signA;
              for (int sB = 0; sB < llistB[OV]; sB++) {
                int addrB1 = abs(listB[OV][k*nmo+l][sB]) - 1;
                int signB = 2 * (listB[OV][k*nmo+l][sB] > 0) - 1;
                int addrB2 = abs(listB[OV][l*nmo+k][sB]) - 1;
                double V_AB = V_A * signB;
                y[addrA1][addrB1] += x[addrA2][addrB2] * V_AB;
                y[addrA1][addrB2] += x[addrA2][addrB1] * V_AB;
                y[addrA2][addrB1] += x[addrA1][addrB2] * V_AB;
                y[addrA2][addrB2] += x[addrA1][addrB1] * V_AB;
              }
            }
          }
        }
      }
    }
    t_B.finish();
    //t_B.print();
    // 2e same spin 3-index part
    Timer t_3A("2e 3-index A");
    t_3A.start();
    for (int j = 0; j < nmo; j++) {
      for (int k = 0; k < j; k++) {
        for (int i = 0; i < nmo; i++) {
          if (i == j || i == k) {
            continue;
          }
          double V =  tei[ind1(ind1(i,i),ind1(j,k))] - tei[ind1(ind1(i,k),ind1(j,i))];
          #pragma omp for schedule(guided, 1)
          for (int sA = 0; sA < llistA[OOV]; sA++) {
            int addrA1 = abs(listA[OOV][(i*nmo+j)*nmo+k][sA]) - 1; // OOV(i,j,k;s)
            int sign = 2 * (listA[OOV][(i*nmo+j)*nmo+k][sA] > 0) - 1;
            int addrA2 = abs(listA[OOV][(i*nmo+k)*nmo+j][sA]) - 1;  // OOV(i,k,j;s)
            double V_signed = V * sign;
            for (int addrB = 0; addrB < nB; addrB++) {
               y[addrA1][addrB] += x[addrA2][addrB] * V_signed;
               y[addrA2][addrB] += x[addrA1][addrB] * V_signed;
            }
          }
        }
      }
    }
    t_3A.finish();
    //t_3A.print();
    Timer t_3B("2e 3-index B");
    t_3B.start();
    #pragma omp for schedule(guided, 1)    
    for (int addrA = 0; addrA < nA; addrA++) {
      for (int j = 0; j < nmo; j++) {
        for (int k = 0; k < j; k++) {
          for (int i = 0; i < nmo; i++) {
            if (i == j || i == k) {
              continue;
            }
            double V =  tei[ind1(ind1(i,i),ind1(j,k))] - tei[ind1(ind1(i,k),ind1(j,i))];
            for (int sB = 0; sB < llistB[OOV]; sB++) {
              int addrB1 = abs(listB[OOV][(i*nmo+j)*nmo+k][sB]) - 1; // OOV(i,j,k;s)
              int sign = 2 * (listB[OOV][(i*nmo+j)*nmo+k][sB] > 0) - 1;
              int addrB2 = abs(listB[OOV][(i*nmo+k)*nmo+j][sB]) - 1;  // OOV(i,k,j;s)
              double V_signed = V * sign;
              y[addrA][addrB1] += x[addrA][addrB2] * V_signed;
              y[addrA][addrB2] += x[addrA][addrB1] * V_signed;
            }
          }
        }
      }
    }
    t_3B.finish();
    //t_3B.print();
    // 2e same spin 4-index part
    Timer t_4A("2e 4-index A");
    t_4A.start();
    for (int i = 0; i < nmo; i++) {
      for (int j = 0; j < i; j++) {
        for (int k = 0; k < i; k++) {
          if (k == j) {
            continue;
          }
          for (int l = 0; l < k; l++) {
            if (l == j) {
              continue;
            }
            double V = tei[ind1(ind1(i,k),ind1(j,l))] - tei[ind1(ind1(i,l),ind1(j,k))]; // (ik||jl)-(il||jk)
            # pragma omp for schedule(guided, 1)
            for (int sA = 0; sA < llistA[OOVV]; sA++) {
              int addrA1 = abs(listA[OOVV][ind2(i,j)*ind2(nmo,0)+ind2(k,l)][sA]) - 1;
              int sign = 2 * (listA[OOVV][ind2(i,j)*ind2(nmo,0)+ind2(k,l)][sA] > 0) - 1;
              int addrA2 = abs(listA[OOVV][ind2(k,l)*ind2(nmo,0)+ind2(i,j)][sA]) - 1;
              double V_signed = V * sign;
              for (int addrB = 0; addrB < nB; addrB++) {
                y[addrA1][addrB] += x[addrA2][addrB] * V_signed;
                y[addrA2][addrB] += x[addrA1][addrB] * V_signed;
              }
            }
          }
        }
      }
    }
    t_4A.finish();
    //t_4A.print();
    Timer t_4B("2e 4-index B");
    t_4B.start();
    # pragma omp for schedule(guided, 1)
    for (int addrA = 0; addrA < nA; addrA++) {
      for (int i = 0; i < nmo; i++) {
        for (int j = 0; j < i; j++) {
          for (int k = 0; k < i; k++) {
            if (k == j) {
              continue;
            }
            for (int l = 0; l < k; l++) {
              if (l == j) {
                continue;
              }
              double V = tei[ind1(ind1(i,k),ind1(j,l))] - tei[ind1(ind1(i,l),ind1(j,k))]; // (ik||jl)-(il||jk)
              for (int sB = 0; sB < llistB[OOVV]; sB++) {
                int addrB1 = abs(listB[OOVV][ind2(i,j)*ind2(nmo,0)+ind2(k,l)][sB]) - 1;
                int sign = 2 * (listB[OOVV][ind2(i,j)*ind2(nmo,0)+ind2(k,l)][sB] > 0) - 1;
                int addrB2 = abs(listB[OOVV][ind2(k,l)*ind2(nmo,0)+ind2(i,j)][sB]) - 1;
                double V_signed = V * sign;
                y[addrA][addrB1] += x[addrA][addrB2] * V_signed;
                y[addrA][addrB2] += x[addrA][addrB1] * V_signed;
              }
            }
          }
        }
      }
    }
    t_4B.finish();
    //t_4B.print();
    }
  } else {
  }
  return 0;
}
