%% EX_3_MainDriver_real_true.m
% 
% This script sets up the Respirometry problem from on/off data
% We use computational tools as described in:
%
%   Cho, Chung, and Pendar. "Computational Tools for Inversion and 
%   Uncertainty Estimation in Respirometry". 2021
%
% Run inverse problems: 
%     inverse problems to reconstruct x(t) from experimental h(t)
%
% Here we compare the reconstructions with true 
%
% T.Cho, J.Chung, and H.Pendar, 03/31/21

%% -------- Set x and h --------
% size of signals is 13413
obsize = 13413;
params = Ex3_GetSetup_Real_true(obsize);

%% -------- Plot figures of h(impulse) and b(observation) --------
% Figure 11 in paper

figure,
subplot(1,3,[1 2])
plot(params.b,'linewidth', 1.5),
legend('CO_2 Observation','location','northeast')
axis([1 params.n -0.1 1.2*max(params.b(:))]);

subplot(1,3,3)
hold on, 
plot(params.delay+1:params.delay+params.supp,...
    params.ht(params.delay+1:params.delay+params.supp),...
    'r', 'linewidth', 1.5),
hold off
xlim([params.delay+1, params.delay+params.supp])
legend('Support of h')
set(gcf,'Position',[300 300 600 200])

%% -------- Test linear systems w/ experimental h --------
% Figure 12 in paper
Lin = Ex3_linear_Real(params);