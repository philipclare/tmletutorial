# TMLE Tutorial
R code and data for TMLE tutorial.

This repository contains data and a number of snippets of R code, used in the TMLE tutorial currently under review at Statistics in Medicine.

The repository includes a number of files of R-code:
* data.creation.R - Code for creating the longitudinal dataset used in example analyses. Data creation is done using the package 'simcausal' (1).
* cross-sectional.R - Code for cross-sectional TMLE analysis, both manually and using the package 'tmle' (2).
* long-single-y.R - Code for longitudinal TMLE with a single outcome measurement, both manually and using the package 'ltmle' (3).
* long-repated-y.R - Code for longitudinal TMLE with a repeated outcome measurement, both manually and using the package 'ltmle' (3).

The longitudinal dataset, ldata.RData, is also included in the repository.

1. Sofrygin O, van der Laan Mark J, Neugebauer R. simcausal R Package: Conducting Transparent and Reproducible Simulation Studies of Causal Effect Estimation with Complex Longitudinal Data. Journal of Statistical Software. 2017;81(2):1-47.
2. Gruber S, van der Laan MJ. tmle: An R Package for Targeted Maximum Likelihood Estimation. 2012. 2012;51(13):35.
3. Lendle SD, Schwab J, Petersen ML, van der Laan MJ. ltmle: An R Package Implementing Targeted Minimum Loss-Based Estimation for Longitudinal Data. Journal of Statistical Software. 2017;81(1):21.
