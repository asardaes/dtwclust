This folder contains all the scripts to execute a set of timing experiments for `dtwclust`.

The `CreateCharTraj.m` script should be run with MATLAB or Octave.
It will download the data and convert it to a series of CSV files that can be imported into `R`.

The `00-main.R` file is the entry point for the experiments.
It loads the necessary packages and runs the other scripts.

The experiments are pretty much independent of each other,
they only need the necessary packages and the data to be loaded.
The data is imported into `R` by `10-read-csv.R`.

The scripts themselves contain some more information and comments.
