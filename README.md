# RBM_TRBM
RBM, TRBM and neural metrics

Models are presented in the article:  [Blindfold learning of an accurate neural metric, Gardella, Marre and Mora, 2017, arXiv:1710.05067](https://arxiv.org/abs/1710.05067#).

This repository allows you to learn the Restricted Boltzmann Machine (RBM) and Temporal Restricted Boltzmann Machine (TRBM).
It also explains how to compute neural metrics based on the RBM and TRBM.
Examples are provided in scripts `MAIN_RBM.m`, `MAIN_TRBM.m` and `MAIN_metrics.m`.

## Warning
The code uses .mex functions, allowing Matlab to run C code. The .mex files must be compiled before using them on
 a new computer, by going into folder `Mex_functions` and running script `COMPILE_mex_files.m`.
