function Lin = Ex2_linear_Real(output)
%% Linear solvers
% Tikhonov(backslash)
Lin.opts_bs.lambda = 0.1;
Lin.opts_bs.L = output.L;
[Lin.X_bs,Lin.INFO_bs] = linsolve_respirometry(output.Hfull,output.b,'\',Lin.opts_bs);
% HyBR
Lin.opts_hybr.mu = sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_HyBR,Lin.INFO_HyBR] = linsolve_respirometry(output.H,output.b,'HyBR',Lin.opts_hybr);
% FLSQR
Lin.opts_flsqr.mu = sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_FLSQR,Lin.INFO_FLSQR] = linsolve_respirometry(output.H,output.b,'FLSQR',Lin.opts_flsqr);
% FISTA
Lin.opts_fista.lambda = 0.002;
Lin.opts_fista.mu = sum(output.b)/output.n *ones(output.n,1); 
Lin.opts_fista.mu = 0*sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_FISTA,Lin.INFO_FISTA] = linsolve_respirometry(output.H,output.b,'FISTA',Lin.opts_fista);

min_fig = min( [min(Lin.X_bs), min(Lin.X_HyBR), min(Lin.X_FLSQR), min(Lin.X_FISTA), min(output.b)] );
max_fig = max( [max(Lin.X_bs), max(Lin.X_HyBR), max(Lin.X_FLSQR), max(Lin.X_FISTA), max(output.b)] );
% Figures
figure, 
hold on
    plot(output.b);
    plot(Lin.X_bs); 
    plot(Lin.X_HyBR,'--k');
    plot(Lin.X_FLSQR); 
    plot(Lin.X_FISTA);
    h1=fill([1 1 output.delay output.delay], [min_fig max_fig max_fig min_fig] , 'r',...
        [output.n-output.delay+1 output.n-output.delay+1 output.n output.n],...
        [min_fig max_fig max_fig min_fig] , 'r');
    set(h1,'FaceAlpha',0.1);
    set(h1,'Edgecolor','none');
    h2=fill([output.delay+1 output.delay+1 output.delay+output.supp output.delay+output.supp], ...
        [min_fig max_fig max_fig min_fig] , 'b',...
        [output.n-output.supp-output.delay+1 output.n-output.supp-output.delay+1 output.n-output.delay+1 output.n-output.delay+1], ...
        [min_fig max_fig max_fig min_fig] , 'b');
    set(h2,'FaceAlpha',0.1);
    set(h2,'Edgecolor','none');
legend('observation','Tikhonov','HyBR','FLSQR','FISTA','delay','supp','location','northwest')
title(['Estimations - ',output.bndry])
% ylim([0, 1]);
set(gcf,'Position',[300 300 1200 400])
hold off


