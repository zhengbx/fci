CC=g++
CFLAGS=-c -Wall -funroll-loops
DEBUG=-g
OMP=-fopenmp
SRC=$(wildcard *.cpp)
OBJ=$(SRC:%.cpp=%.o)
OPT=-O3 -std=c++0x
#OPT=-O0
MKLROOT=/opt/intel/mkl
MKLL=-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
MKLC=-m64 -I$(MKLROOT)/include

EXE=fci

all:$(EXE)

$(EXE):$(OBJ) makefile
	$(CC) $(OBJ) $(DEBUG) $(OMP) $(LIBS) $(MKLL) -o $@
%.o:%.cpp %.h makefile
	$(CC) $< $(DEBUG) $(CFLAGS) $(OMP) $(MKLC) $(OPT) 
clean:
	-rm -f $(EXE) $(OBJ)
