#!/bin/bash

if [[ ! -d "data" ]]; then
    mkdir data
fi

make main

if [[ $? ]]; then
    ./simulation.out 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.005 0 0 0 "data/run1.txt"
    ./simulation.out 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.005 0 10000000000 2 "data/run2.txt"
    ./simulation.out 0.005 0.000025 0.008 0.25 0.05 200 0.5 0.005 0 0 0 "data/run3.txt"
else
    echo "Compilation failed."
fi