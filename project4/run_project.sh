#!/bin/bash

main="run_ising.out"

if [[ ! -d "data" ]]; then
    mkdir data
fi

if [[ (! -d "data/test")  ||  (! -d "data/burn_in")  || (! -d "data/approx_distr")  ||  (! -d "data/pahse") ]]; then
    mkdir data/test data/burn_in data/approx_distr data/phase
fi

make all

./$main

./plot.sh