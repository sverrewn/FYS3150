#!/bin/bash

cmd="save"

if [[ $# -gt 0 ]]; then
    option="$1"
    case $option in
        -d|--delete)
        cmd="del"
        ;;
        -s|--save)
        cmd="save"
        ;;
        -h|--help|*)
        echo "./run.sh [{-d | -s}]"
        echo "-d, --delete  delete all .dat files after running script"
        echo "-s, --save    save all .dat files in the data folder"
        echo "Default behaviour is is -s"
        exit 0
        ;;
    esac
fi

F1="poisson_exact.out"
F2="general_tridiag.out"

dir1="plots"
dir2="data"

if [ ! -d "$dir1" ]; then
    mkdir $dir1
fi

if [[ ( ! -f "$F1" ) || ( ! -f "$F2" ) ]]; then
    make release
fi

nums=(10 100 1000 10000)

for val in ${nums[@]}; do
    ./$F1 $val
    ./$F2 $val
done

python plot_pr2.py poisson_exact_1000.dat
python plot_pr7.py $(ls general_tridiag_*.dat)
python plot_pr8.py $(ls poisson_exact_*.dat) $(ls general_tridiag_*.dat)

case $cmd in
    "del")
    rm *.dat
    ;;
    "save")
    if [[ ! -d "$dir2" ]]; then
        mkdir data
    fi
    mv *.dat data
    ;;
esac
