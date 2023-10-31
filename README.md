# modeling-conditionals-empirical

This repository contains the code for Chapter 7 of my Ph.D thesis (see [dissertation](https://osnadocs.ub.uni-osnabrueck.de/handle/ds-202310249948)), in which
I apply our RSA-model for the communication with conditionals (see [paper](https://semprag.org/index.php/sp/article/view/sp.15.13/3055)) on the empirical data from a behavioral experiment that we run (see [this repository](https://github.com/brittaGrusdt/communicating-uncertain-beliefs-conditionals)).



- R/single-model-run.R
Run the RSA-model once with a single set of fixed parameters.


- R/fit-rsa-model.R
Fit the RSA model to the empirical data. On top of this file, specify speaker type, parameters to be fitted etc. , all configured in the config.yml file.

- R/fit-data-independent-conditions.R and R/fit-data-dependent-conditions.R
ZOIB-models fitted to empirical data from the PE-task of the experiment


- R/independent-pe-task-dirichlet.R and R/dependent-pe-task-dirichlet.R
Dirichlet-models fitted to empirical data from the PE-task of the experiment


- R/plots-behavioral-data.R
Plots of behavioral data as they appear in my Ph.D thesis
