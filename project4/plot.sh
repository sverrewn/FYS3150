#!/bin/bash

if [[ ! -d "figs" ]]; then
    mkdir figs
fi

python scripts/plot4.py $(ls data/test)
python scripts/plot5.py $(ls data/burn_in)
python scripts/plot6.py $(ls data/approx_distr)

python scripts/plot8.py $(ls data/phase/L40*) $(ls data/phase/L60*) $(ls data/phase/L80*) $(ls data/phase/L100*)
