function Lin = Ex3_linear_Real(output)
%% Linear solvers

% HyBR
Lin.opts_hybr.mu = sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_HyBR,Lin.INFO_HyBR] = linsolve_respirometry(output.H,output.b,'HyBR',Lin.opts_hybr);
% FLSQR
Lin.opts_flsqr.mu = sum(output.b)/output.n *ones(output.n,1); 
Lin.opts_flsqr.mu = 0*sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_FLSQR,Lin.INFO_FLSQR] = linsolve_respirometry(output.H,output.b,'FLSQR',Lin.opts_flsqr);
% FISTA
Lin.opts_fista.lambda = 0.002;
Lin.opts_fista.mu = sum(output.b)/output.n *ones(output.n,1); 
Lin.opts_fista.mu = 0*sum(output.b)/output.n *ones(output.n,1); 
[Lin.X_FISTA,Lin.INFO_FISTA] = linsolve_respirometry(output.H,output.b,'FISTA',Lin.opts_fista);

min_fig = min( [min(Lin.X_HyBR(1001:5000)), min(Lin.X_FLSQR(1001:5000)), min(Lin.X_FISTA(1001:5000)), min(output.b(1001:5000))] );
max_fig = max( [max(Lin.X_HyBR(1001:5000)), max(Lin.X_FLSQR(1001:5000)), max(Lin.X_FISTA(1001:5000)), max(output.b(1001:5000))] );


%% Plot figures
zoom1 = 3301;
zoom2 = 5063;
figure,
subplot(2,1,1)
hold on
    plot(output.b);
    plot(output.x_exact,'color',[0 0 0] +0.5,'linewidth',1.2);
    plot(Lin.X_HyBR,'color',[0.8500 0.3250 0.0980]);
    plot(Lin.X_FLSQR,'color',[0.9290 0.6940 0.1250]); 
    plot(Lin.X_FISTA,'color',[0.4660 0.6740 0.1880]);
    plot(zoom1:1:zoom2, max_fig*ones(1,zoom2-zoom1+1),'-k','linewidth',3);
    plot(zoom1:1:zoom2, min_fig*ones(1,zoom2-zoom1+1),'-k','linewidth',3);
    plot(zoom1*ones(1,2), [min_fig max_fig],'-k','linewidth',3);
    plot(zoom2*ones(1,2), [min_fig max_fig],'-k','linewidth',3);
    set(gca,'fontsize',20)
    title('CO_2 Reconstructions (Full size)','fontsize',20)
    ylim([min_fig, max_fig]);
    xlim([1 output.n]);
hold off

subplot(2,1,2)
hold on,
    plot(output.b,'linewidth',1.5);
    plot(output.x_exact,'color',[0 0 0] +0.5,'linewidth',1.8);
    plot(Lin.X_HyBR,'color',[0.8500 0.3250 0.0980],'linewidth',1.5);
    plot(Lin.X_FLSQR,'color',[0.9290 0.6940 0.1250],'linewidth',1.5);
    plot(Lin.X_FISTA,'color',[0.4660 0.6740 0.1880],'linewidth',1.5);
    set(gca,'fontsize',20)
    legend('Observation','True','HyBR','FLSQR-R','FISTA','location','northeast','fontsize',15)
    title('CO_2 Reconstructions (Zoom in)','fontsize',20)
    ylim([min_fig, max_fig]);
    xlim([zoom1 zoom2]);
hold off

arh1=annotation('arrow', [0.32,0.13 ],[0.58,0.46] );
arh2=annotation('arrow', [0.425,0.905 ],[0.58,0.46] );
arh1.LineWidth = 3; arh2.LineWidth = arh1.LineWidth;
arh1.Color = [0,0,0]; arh2.Color = arh1.Color;

set(gcf,'Position',[300 300 1200 800])

