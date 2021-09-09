#!/bin/bash

F1="general_tridiag.out"

if [ ! -f "$F1" ]; then
    make release
fi

nums=(10 100 1000 10000 100000)

for val in ${nums[@]}; do
    ./$F1 $val
done