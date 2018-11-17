#!/bin/bash

# Run Lorentz 96 and compute covariance with 1st and 2nd order derivatives

# ARG1: dimension
# ARG2: forcing
# ARG3: timesteps


time env OMP_NUM_THREADS=1 ./plorenz --tensor1 6 10 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor2 6 10 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor3 6 10 40000

./plot.sh

cp *.png figures/6_10

time env OMP_NUM_THREADS=1 ./plorenz --tensor1 6 3.6 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor2 6 3.6 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor3 6 3.6 40000

./plot.sh

cp *.png figures/6_3_6

time env OMP_NUM_THREADS=1 ./plorenz --tensor1 7 2.0 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor2 7 2.0 40000
time env OMP_NUM_THREADS=1 ./plorenz --tensor3 7 2.0 40000

./plot.sh

cp *.png figures/7_2_0

time ./plorenz --tensor1 7 4.4 40000
time ./plorenz --tensor2 7 4.4 40000
time ./plorenz --tensor3 7 4.4 40000

./plot.sh

cp *.png figures/7_4_4
