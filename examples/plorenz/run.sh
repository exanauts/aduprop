#!/bin/bash

# Run Lorentz 96 and compute covariance with 1st and 2nd order derivatives

# ARG1: dimension
# ARG2: forcing
# ARG3: timesteps

./plorenz --tensor1 $1 $2 $3 | grep 'COV: ' | sed 's/[COV:|\[|,]//g' | sed 's/\]//' | awk '{$1=""; print $0}' > data_cov1

./plorenz --tensor2 $1 $2 $3 | grep 'COV: ' | sed 's/[COV:|\[|,]//g' | sed 's/\]//' | awk '{$1=""; print $0}' > data_cov2

./plorenz --tensor3 $1 $2 $3 | grep 'COV: ' | sed 's/[COV:|\[|,]//g' | sed 's/\]//' | awk '{$1=""; print $0}' > data_cov3

./plorenz $1 $2 $3 | grep 'X: ' | sed 's/[X:|\[|,]//g' | sed 's/\]//' | awk '{$1=""; print $0}' > datas
