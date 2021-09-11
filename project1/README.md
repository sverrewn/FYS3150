# Project 1

The easiest way to compile and test the project is to run the script run.sh`./run.sh [-s | -d]`. The script will by default save (the `-s` flag) all data files generated in the sub folder data/ and all plots in the sub folder plots/. If you do not want to save the data files or the plots, the flag `-d` (or `--delete`) will remove these after running the code. The script might need permission to run, which is easily solved with `chmod a+x run.sh`



To manually run the project, compile with `make release` which will build the files with the `-Ofast` compiler flag. You can also run `make test` which compiles without optimisations  and with `-Wall -Wextra -Wpedantic`

