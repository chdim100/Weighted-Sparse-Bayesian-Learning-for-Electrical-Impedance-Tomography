# Weighted-Sparse-Bayesian-Learning-for-Electrical-Impedance-Tomography
A MATLAB code package to perform Weighted Sparse Bayesian Learning for EIT using Bound-Optimization. Allows parametrization, compares quantitatively and qualitatively with Regularization (L2 and L1) and EM-based Block Sparse Bayesian Learning

Requires the EIDORS package libary (version 3.9 or later) to be executed. **https://eidors3d.sourceforge.net/download.shtml** 
Uses experimental data from https://zenodo.org/record/1203914#.Y5MGbXZByUk and in-vivo thoracic data from EIDORS. 

Supports only 2D EIT imaging, with 3D imaging planned as a future capability. 

Run the Main.m file, with proper selection of the cases, perturbations and measurement patterns (see comments in file)

A Python demo is also included for Weighted Sparse Bayesian Learning which will be regularly updated. 

Please cite to:  
Dimas, Christos, Vassilis Alimisis, and Paul P. Sotiriadis, "Electrical Impedance Tomography using a Weighted Bound-Optimization Block Sparse Bayesian Learning Approach" 2022 IEEE 22th International Conference on Bioinformatics and Bioengineering (BIBE). IEEE, 2022 DOI: 10.1109/BIBE55377.2022.00058 
