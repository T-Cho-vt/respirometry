%% startup_Respirometry.m
%
% startup file for the respirometry codes

% Required packages
addpath(genpath('RestoreTools'))
addpath(genpath('IRtools-master'))
addpath(genpath('FISTA-master'))
addpath(genpath('SIAMBookCodes-master'))
addpath(genpath('genHyBR-master'))

% Classes, functions, and data for respirometry experiments
addpath(genpath('Class_and_Function'))
addpath(genpath('Real_data'))

% Scripts for MainDrivers
addpath(genpath('Ex1_Simulation'))
addpath(genpath('Ex2_real'))
addpath(genpath('Ex3_real_true'))