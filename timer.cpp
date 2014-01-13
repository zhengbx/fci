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

#include "timer.h"
#include <iostream>
#include <cstdio>

int Timer::print() {
  printf("Time for %-20s%10.2f s\n", (name + ":").c_str(), data);
  return 0;
}

double Timer::time() {
  return data;
}

int Timer::start() {
  gettimeofday(&t_start, NULL);
  return 0;
}

int Timer::finish() {
  timeval t_end;
  gettimeofday(&t_end, NULL);
  data += (double)(t_end.tv_sec-t_start.tv_sec)
      +((double)(t_end.tv_usec-t_start.tv_usec))/1000000.;
  return 0;
}

