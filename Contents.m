%% Contents.m
%   This folder contains MATLAB code to accompany the paper:
%
%   Cho, Chung, and Pendar. "Computational Tools for Inversion and 
%   Uncertainty Estimation in Respirometry". 2021
%
%   The MainDriver's codes require the following packages:
%
%       RestoreTools by James Nagy, Katrina Palmer, and Lisa Perrone
%             http://www.mathcs.emory.edu/~nagy/RestoreTools/index.html
%
%       IRTools by James Nagy, Sivia Gazzola, and Per Christian Hansen
%             https://github.com/silviagazzola/IRtools
%
%       genHyBR by Julianne Chung and Arvind Saibaba
%             https://github.com/juliannechung/genHyBR
%
%       Codes from the SIAM Book by Johnathan M. Bardsley
%       "Computational Uncertainty Quantification for Inverse Problems"
%             https://github.com/bardsleyj/SIAMBookCodes  
%
%       FISTA, a MATLAB implementation by Tiep Vu
%             https://github.com/tiepvupsu/FISTA
%
%   First run startup_Respirometry.m for setting paths to the above packages.
%
%   Cho, Chung, and Pendar (2021)

%% MainDrivers for each numerical experiments of inverse problems
%
%   EX_1_MainDriver_sim.m           Sets up and runs a linear inverse
%                                   problem corresponding to Sections  
%                                   "Linear Respirometry Reconstruction"
%                                   and "Nonlinear Respirometry
%                                   Reconstruction" with simulations.
%                                   Generates Figures 1-10
% 
%   EX_2_MainDriver_real.m          Sets up and runs real blind
%                                   deconvolution corresponding to
%                                   "Nonlinear case study: Abdominal
%                                   pumping and CO2 emission in a darkling
%                                   beetle" and generating Figures 13, 14
% 
%   Ex_3_MainDriver_real_true.m     Sets up and runs real linear case study
%                                   generating Figures 11, 12
