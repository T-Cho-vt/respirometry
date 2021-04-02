%% EX_2_MainDriver_real.m
% 
% This script sets up the Respirometry problem for real data
% We use computational tools as described in:
%
%   Cho, Chung, and Pendar. "Computational Tools for Inversion and 
%   Uncertainty Estimation in Respirometry". 2021
%
% Run inverse problems: 
%     nonlinear inverse problems to reconstruct h(t) and x(t)
%     where we assume that the support and delay of h(t) are known
%
% T.Cho, J.Chung, and H.Pendar, 03/31/21

%% -------- Set x and h --------
% Actual size of observation from experiment is 12000 
% from sampling with 10Hz

obsize = 12000;
params = Ex2_GetSetup_Real(obsize);

%% -------- Test Blind Deconvolution --------
% Parameter setup for blind deconvolution with real observations
Params = Ex2_Param_Setup_Real(params);
% Note that we assume delay of h is 213 (e.g., 21.3 seconds)
Params.delay = 213;

% Test 3 different initial h0

% uniform 
h0 = ones(params.supp,1)/params.supp; % initial guess
[blind_x_uni,blind_h_uni,blind_Hs_uni,blind_res_uni,blind_flag_uni] = Ex2_blind_Real(h0,Params);

% gamma dist initial guess
[~,h0,~]=h_Gamma([3 3], 25, 10, 0.0001); h0= h0/sum(h0); h0 = h0(1:params.supp);
[blind_x_gamma,blind_h_gamma,blind_Hs_gamma,blind_res_gamma,blind_flag_gamma] = Ex2_blind_Real(h0,Params);

[~,h0,~]=h_Gamma([1 .5], 25, 10, 0.0001); h0= h0/sum(h0); h0 = h0(1:params.supp);
[blind_x_gamma2,blind_h_gamma2,blind_Hs_gamma2,blind_res_gamma2,blind_flag_gamma2] = Ex2_blind_Real(h0,Params);

% save results
% save('Ex2_bd_full_results.mat','Params',...
%     'blind_x_uni','blind_h_uni','blind_Hs_uni','blind_res_uni','blind_flag_uni',...
%     'blind_x_gamma','blind_h_gamma','blind_Hs_gamma','blind_res_gamma','blind_flag_gamma',...
%     'blind_x_gamma2','blind_h_gamma2','blind_Hs_gamma2','blind_res_gamma2','blind_flag_gamma2')

% load results
% load('Ex2_bd_full_results.mat')

%% Summary
% Figure 13 in paper

% We compare results for first 400 seconds
sig_size = 4000;

% Load abdominal pumping observations
Abdo = textread('AP_30Hz.txt');
% Resample from 30Hz to 10Hz
Abdo = resample(Abdo(:,2),1,3);
Abdo = Abdo(1:sig_size);

% Load photo of living organism in respirometry chamber
Insect_color = imread('Insect_color.jpg');
Insect_color = Insect_color(361:800,301:end-360,:);
Insect_red = imread('Insect_red.jpg');
Insect_red = Insect_red(24:end-51,81:end-80,:);
Insect_red = Insect_red(210:465,175:end-209,:);

recon_sig = blind_x_gamma(1:sig_size);
observ_sig = params.b(1:sig_size);

figure, 
    
    subplot(6,2,1), hold on
        imshow(Insect_color),
        text(0,-70,'(a)','fontsize',20),
    hold off
    
    subplot(6,2,2), hold on
        imshow(Insect_red),
        text(0,-30,'(b)','fontsize',20),
    hold off
    
    subplot(6,2,[3 4]), hold on,
        plot(Abdo,'color',[0.8500, 0.3250, 0.0980], 'linewidth', 1.5),
        text(0,.95,'(c)','fontsize',20),
    hold off
    legend('Abdominal Movement','location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis([1 sig_size 0.1 1.2*max(Abdo(:))]);
    
    subplot(6,2,[5 6]), hold on,
        plot(observ_sig,'color',[0,0,0]+0.45,'linewidth', 1.5),
        text(0,0.8,'(d)','fontsize',20),
    hold off
    legend('Recorded CO$_2$','location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis([1 sig_size -0.01 1.2*max(observ_sig(:))]);
    
    subplot(6,2,[7 8]), hold on,
        plot(recon_sig,'color',[0, 0.4470, 0.7410],'linewidth', 1.2),
        text(0,1.75,'(e)','fontsize',20),
    hold off
    legend('Corrected CO$_2$','location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    axis([1 sig_size -0.01 1.2*max(recon_sig(:))]);
    
    subplot(6,2,[9 10 11 12]), hold on,
        plot(observ_sig,'color',[0,0,0]+0.45,'linewidth', 1.5),
        plot(Abdo+0.85,'color',[0.8500, 0.3250, 0.0980],'linewidth', 1.5),
        plot(recon_sig,'color',[0, 0.4470, 0.7410],'linewidth', 1.2),
        text(0,1.65,'(f)','fontsize',20),
    hold off
    set(gca,'ytick',[])
    set(gca,'xtick',[0 sig_size/2 sig_size],'xticklabel',{'0','200','400(s)'},'fontsize',15)
    
set(gcf,'Position',[50 100 800 1000])


%% Plot reconstructions of h(t)
% Figure 14 in paper
h_unif = ones(params.supp,1)/params.supp; % initial guess
n2 = floor(params.supp/2);
h0 = linspace(2/params.supp,0,n2)'; h0 = [h0; zeros(params.supp-n2,1)];% initial guess
h_desc = fliplr(h0);
h_exp = params.ht(params.delay+1:params.delay+params.supp)/sum(params.ht); % initial guess
[~,h_gam1,~]=h_Gamma([3 3], 25, 10, 0.0001); h_gam1= h_gam1/sum(h_gam1); h_gam1 = h_gam1(1:params.supp);
[~,h_gam2,~]=h_Gamma([1 .5], 25, 10, 0.0001); h_gam2= h_gam2/sum(h_gam2); h_gam2 = h_gam2(1:params.supp); figure, plot(h_gam2);
subplot(1,2,1)
hold on
    plot(h_unif,'linewidth',2.5,'color',[0.9290, 0.6940, 0.1250])
    plot(h_gam1,'--','linewidth',2.5,'color',[0.8500, 0.3250, 0.0980])
    plot(h_gam2,':','linewidth',2.5,'color',[0.4940, 0.1840, 0.5560])
hold off
legend('Flat line','Gamma1','Gamma2','fontsize',20)
title('Initital h','fontsize',25)
xlim([1 params.supp])
ylim([0 0.07])
set(gca,'fontsize',20)

subplot(1,2,2)
hold on
    plot(blind_h_uni,'linewidth',2.5,'color',[0.9290, 0.6940, 0.1250])
    plot(blind_h_gamma,'--','linewidth',2.5,'color',[0.8500, 0.3250, 0.0980])
    plot(blind_h_gamma2,':','linewidth',2.5,'color',[0.4940, 0.1840, 0.5560])
hold off
xlim([1 params.supp])
ylim([0 0.07])
legend('Flat line - est','Gamma1 - est','Gamma2 - est','fontsize',20)
title('Reconstructed h','fontsize',25)
set(gca,'fontsize',20)

set(gcf,'Position',[300 300 1200 400])
