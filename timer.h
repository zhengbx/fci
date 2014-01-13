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

#ifndef TIME
#define TIME

#include "sys/time.h"
#include <string>

class Timer {
  double data;
  std::string name;
  timeval t_start;
public:
  Timer(std::string _name = ""): data(0), name(_name) {};
  int print();
  int start();
  int finish();
  double time();
};

#endif
