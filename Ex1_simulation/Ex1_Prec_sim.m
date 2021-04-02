function output_prec = Ex1_Prec_sim(n,output)
%% Clustering for different tau
% Use GCV to select tau
tau1 = GCVforSVD(fft(output.Hfull(:,1)), fft(output.b)); 
output_prec.tau1 = tau1;
tau2 = tau1*10;
output_prec.tau2 = tau2;
tau3 = tau1/10;
output_prec.tau3 = tau3;

M1 = precMatrix(output.Hfull(:,1),tau1); M1full = []; 
M2 = precMatrix(output.Hfull(:,1),tau2); M2full = [];
M3 = precMatrix(output.Hfull(:,1),tau3); M3full = [];

I = eye(n);
for i = 1:n
    M1full(:,i) = M1*I(:,i); 
    M2full(:,i) = M2*I(:,i); 
    M3full(:,i) = M3*I(:,i);
end
invM1H = M1full \ output.Hfull;
invM2H = M2full \ output.Hfull;
invM3H = M3full \ output.Hfull;

% figure, 
% subplot(1,3,1), imagesc(output.Hfull)
% subplot(1,3,2), imagesc(M1full)
% subplot(1,3,3), imagesc(invM1H)

[U_unc,s_unc,V_unc] = svd(output.Hfull);
[U_prec1,s_prec1,V_prec1] = svd(invM1H);
[U_prec2,s_prec2,V_prec2] = svd(invM2H);
[U_prec3,s_prec3,V_prec3] = svd(invM3H);

figure,semilogy(diag(s_unc),'b*'),hold on, 
semilogy(diag(s_prec3),'--k','linewidth',2)
semilogy(diag(s_prec1),'--r','linewidth',2)
semilogy(diag(s_prec2),'--c','linewidth',2)


hold off; axis([10 n 1e-5 8])
legend('Unprecon',...
    sprintf('Precon tau = %.1d',tau3),...
    sprintf('Precon tau = %.1d',tau1),...
    sprintf('Precon tau = %.1d',tau2),...
    'location','southwest')
ylabel('singular values')

set(gcf,'Position',[300 300 600 200])

%% Regularized PCGLS for right and left
% nLevel: 0.01(high), 0.0075(middle), 0.005(low,default)
nlevelvec = [0.01, 0.0075, 0.005];
for i = 1:length(nlevelvec)
    output = Ex1_GetSetup_sim(n,nlevelvec(i));

    Efun = @(alpha) norm( (output.Hfull'*output.Hfull + alpha.^2*output.Q'*output.Q) \ (output.Hfull'*output.b) - output.x_exact );
    alpha_reg = fminbnd(Efun,0,1);
    alpha_reg = sqrt(alpha_reg);
    Hfull_reg = [output.Hfull; alpha_reg^2*output.Q];
    b_reg = [output.b; zeros(n,1)];
    output_prec.alpha_reg = alpha_reg;

    % CLGS
    [x1, k1, Rnrm1, Xnrm1, Enrm1] = CGLS(output.Hfull, output.b, output.x_exact*0, 200, [], output.x_exact);
    [x1_reg, k1_reg, Rnrm1_reg, Xnrm1_reg, Enrm1_reg] = CGLS(Hfull_reg, b_reg, output.x_exact*0, 200, [], output.x_exact);

    % L-CGLS
    [xL1,kL1,RnrmL1,XnrmL1,EnrmL1] = CGLS(M1full\output.Hfull, M1full\output.b, output.x_exact*0, 200,[],output.x_exact);
    [xL1_reg,kL1_reg,RnrmL1_reg,XnrmL1_reg,EnrmL1_reg] = CGLS([M1full\output.Hfull; alpha_reg^2*output.Q], [M1full\output.b; zeros(n,1)], output.x_exact*0, 200,[],output.x_exact);

    % R-CGLS
    [xp1, kp1, Rnrmp1, Xnrmp1, Enrmp1] = PCGLS(output.H, M1, output.b, output.x_exact*0, 200, [], output.x_exact);
    [xp1_reg, kp1_reg, Rnrmp1_reg, Xnrmp1_reg, Enrmp1_reg] = PCGLS(Hfull_reg, M1, b_reg, output.x_exact*0, 200, [], output.x_exact);
    %%
    figure,
    hold on, 
     semilogy(Enrm1,'b','linewidth',1), 
     semilogy(Enrm1_reg,'--b','linewidth',1)
     semilogy(EnrmL1,'k','linewidth',1),
     semilogy(EnrmL1_reg,'--k','linewidth',1),
     semilogy(Enrmp1,'r','linewidth',1), 
     semilogy(Enrmp1_reg,'--r','linewidth',1), 
     plot(k1,Enrm1(k1),'xb','linewidth',2)
     plot(k1_reg,Enrm1_reg(k1_reg),'sb','linewidth',2)
     plot(kL1,EnrmL1(kL1),'xk','linewidth',2)
     plot(kL1_reg,EnrmL1_reg(kL1_reg),'sk','linewidth',2)
     plot(kp1,Enrmp1(kp1),'xr','linewidth',1.5)
     plot(kp1_reg,Enrmp1_reg(kp1_reg),'sr','linewidth',1.5)
    hold off, 
    legend('CGLS','CGLS-Tik','Left-PCGLS','Left-PCGLS-Tik','Right-PCGLS','Right-PCGLS-Tik',...
        'NumColumns',3)
    axis([1 length(Enrm1) 0.05 1])
    title(sprintf('nLevel: %.2f %%',output.nLevel*100))
    set(gcf,'Position',[300 300 600 120])
end