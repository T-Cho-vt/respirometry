%% EX_1_MainDriver_sim.m
%
% This script sets up and run Respirometry problem for simulation data.
% We use computational tools as described in:
%
%   Cho, Chung, and Pendar. "Computational Tools for Inversion and 
%   Uncertainty Estimation in Respirometry". 2021
%
% Run inverse problems:
%  1. linear inverse problems and uncertainty quantifcation
%  2. acceleration iterative methods
%  3. nonlinear inverse problems to reconstruct h(t) and x(t)
%     where we assume that the support and delay of h(t) are known
%
% T.Cho, J.Chung, and H.Pendar, 03/31/21

%% -------- Size of signal and Set noise level --------

% size of signal and impulse response function
% Note that figures from this script are optimized for n=512;
n= 512;

% nLevel: 0.01(high), 0.0075(middle), 0.005(low,default)
nLevel = 0.005;

%% -------- Set x, h, and params --------
%
% params is structure containing all variables and parameters 
% for linear inverse problems

params = Ex1_GetSetup_sim(n,nLevel);

%% -------- Plot reduced system ---------
% Figure 1 in paper. 
Ex1_Fig_reduced_system;

%% -------- Plot figures of x, b, and h --------
% Figure 4 in paper
figure, 
plot(params.x_exact,'linewidth', 1.5)
hold on, plot(params.b,'r','linewidth', 1.5), hold off
legend('true','observed','location','northeast')
set(gcf,'Position',[300 300 600 120])

figure,
plot(params.ht, 'linewidth', 1.5), hold on, 
plot(params.delay:1:params.delay+params.supp,...
    params.ht(params.delay:1:params.delay+params.supp),...
    '--r', 'linewidth', 1.5),
hold off
legend('h','support')
axis([1 n 0 1.2*max(params.ht(:))]);
set(gcf,'Position',[300 300 600 120])

%% -------- Test linear systems --------
% Figure 5 in paper
output_lin = Ex1_linear_sim(n,params);

%% -------- Test UQ --------
% Figure 6 in paper
output_UQ = Ex1_UQ_sim(n,params);

%% -------- Test Preconditioning --------
% Figure 2 and 7 in paper
output_prec = Ex1_Prec_sim(n,params);


%% -------- Test Blind Deconvolution --------
% Sets up parameters for blind deconvolution with simulation data
Params = Ex1_Param_Setup_sim(n,params);

% initial guess
h0 = ones(params.supp,1)/params.supp; 

% Runs blind deconvolution
[blind_x,blind_h,blind_Hs,blind_res,blind_flag] = Ex1_blind_sim(h0,Params);

%% -------- Plot reconstructed x and h --------
% Figure 8 and 9 in paper

% Plot reconstructed h and x
view_iter_h = [1 51 101 size(blind_Hs,2)];
figure, 
for i = 1:length(view_iter_h)
    subplot(1,length(view_iter_h),i), 
    hold on
        plot(params.ht(params.delay+1:params.delay+params.supp),'linewidth',1.5)
        plot(blind_Hs(:,view_iter_h(i)),'-.','linewidth',1.5) 
    hold off
    xlim([1,params.supp])
    
    if i == length(view_iter_h)
        title(['k = ', num2str(view_iter_h(i)),' (stop)'])
        axes('position',[.85 .6 .05 .25])
        box on
            index_max_ht = find(params.ht==max(params.ht))-params.delay;
            index_max_blind = find(blind_Hs(:,end) == max(blind_Hs(:,end)));
            index_min = min(index_max_ht,index_max_blind);
            index_max = max(index_max_ht,index_max_blind);
            indexOfInterest = index_min-3:1:index_max+5;
            
            hold on
                plot(params.ht(params.delay+indexOfInterest),'linewidth',1.5)
                plot(blind_Hs(indexOfInterest,view_iter_h(i)),'-.','linewidth',1.5)
            hold off
            
            max_h = max(params.ht(params.delay+index_max_ht),blind_Hs(index_max_blind,end));
            axis([1 length(indexOfInterest) max_h*.5 max_h*1.2])
        box off
    else
        title(['k = ', num2str(view_iter_h(i)-1)])
    end
    
    if i == 1
       legend('true','h_k') 
    end
end
set(gcf,'Position',[300 300 800 150])

% from final h, construct x
output2 = params;
output2.ht = [zeros(params.delay,1);blind_h;zeros(n-params.supp-params.delay,1)];
output2.Hfull = matrixH(output2.ht,Params.bndry);
output2.H = matrixHfun(output2.ht,Params.bndry,1);

output_lin2 = Ex1_linear_sim(n,output2);

%% -------- Motivations for constraints on x --------
% Figure 3 in paper

% Showing blinde deconvolution can be stuck with 2-norm reg
Ex1_mot_for_const_sim

%% -------- Inexact delay--------
% Figure 10 in paper

% Runs blind deconvolution experiments with inexact delays
% True delay is 9
% Test: 2,9,16

Ex1_InexactDelay_sim