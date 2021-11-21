#!/bin/bash

python scripts/plot4.py $(ls data/test)
python scripts/plot5.py $(ls data/burn_in)
python scripts/plot6.py $(ls data/approx_distr)