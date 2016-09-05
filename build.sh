#!/bin/bash

g++6.1.0 -Wl,-rpath,$HOME/opt/gcc-6.1.0/lib64 \
--std=c++14 \
-Wall -Wno-reorder -Wno-terminate \
-Ofast -flto -march=native -funroll-loops -fopenmp \
-I. \
-I$HOME/opt/nuSQuIDS/inc -I$HOME/opt/SQuIDS/include -Wno-abi -I$HOME/opt/hdf5-1.8.17/include \
`root-config --cflags` \
-L$HOME/opt/nuSQuIDS/lib -lnuSQuIDS -L$HOME/opt/SQuIDS/lib -lSQuIDS -lgsl -lgslcblas -lm -L$HOME/opt/hdf5-1.8.17/lib -lhdf5 -lhdf5_hl \
`root-config --glibs` -lMathMore \
-Wl,--as-needed \
./MyEarthAtm.cpp ./main.cpp -o atmonu_osc