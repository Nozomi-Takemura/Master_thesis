# Project Title

Masther thesis -Tensor decompositions for latent Dirichlet allocation in topic modeling-

## Description

This repository includes master thesis and the codes written in Python and Matlab used to conduct experiments illustrated in the thesis. Synthetic data generation from latent Dirichlet allocation (LDA) can be done with the Python script in "Python/DataGeneration/", while parameter inference(recovery) with collapsed Gibbs sampling(tensor decomposition) can be performed with the functions written in Matlab located in "Matlab/Function/." The experiments for synthetic data in the thesis could be reproduced with "K2_0423_same_sample_ortho_100_runs_maxite5000.m" and with the codes located in "MatlabCode/ErrorAnalysis/". Finally, all codes utilized for topic modeling for NIPS datasets can be found in "MatlabCode/NIPSdata/".

It is recommended that you use "Python/DataGeneration/loop.py" file to generate a synthetic dataset from LDA since it is more flexible regarding the specification of the parameter of an LDA than the other python scripts. 

### Prerequisites

To use the code in this repository, you need to install Matlab (2017b) and Python3. Furthermore, the tensor decomposition package for Matlab, tensor, needs to be downloaded. 





