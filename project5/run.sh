#!/bin/bash

if [[ ! -d "data" ]]; then
    mkdir data
fi

make main

if [[ $? ]]; then
    ./simulation.out 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.05 0   0 0 "data/run1.bin"
    ./simulation.out 0.005 2.5e-5 0.008 0.25 0.05 200 0.5 0.05 0   1e10 2 "data/run2.bin"
    ./simulation.out 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.05 0.2 1e10 2 "data/run3.bin"
    ./simulation.out 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.05 0.2 1e10 1 "data/run4.bin"
    ./simulation.out 0.005 2.5e-5 0.002 0.25 0.05 200 0.5 0.05 0.2 1e10 3 "data/run5.bin"
else
    echo "Compilation failed."
fi