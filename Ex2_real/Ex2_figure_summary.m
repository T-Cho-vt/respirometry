% Ex 2 figures summary
close all 
delay_bd = 186;
sig_size = 4000;

Abdo = textread('AP_30Hz.txt');
Abdo = resample(Abdo(:,2),1,3);
Abdo = Abdo(1:sig_size);

Insect_color = imread('Insect_color.jpg');
Insect_color = Insect_color(361:800,301:end-360,:);
Insect_red = imread('Insect_red.jpg');
Insect_red = Insect_red(24:end-51,81:end-80,:);
Insect_red = Insect_red(210:465,175:end-209,:);

recon_sig = blind_x_uni(delay_bd+1:sig_size+delay_bd);
observ_sig = params.b(1:sig_size);

nr = 6;
nc = 2;

figure, 
    
    subplot(nr,nc,1), hold on
        imshow(flipud(Insect_color)),
        text(0,500,'(a)','fontsize',20),
    hold off
    
    subplot(nr,nc,2), hold on
        imshow(flipud(Insect_red)),
        text(0,290,'(b)','fontsize',20),
    hold off
    
    subplot(nr,nc,[3 4]), hold on,
        plot(Abdo,'color',[0.8500, 0.3250, 0.0980], 'linewidth', 1.5),
        text(0,.95,'(c)','fontsize',20),
    hold off
    legend('Abdominal Movement','location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    axis([1 sig_size 0.1 1.2*max(Abdo(:))]);
    
    subplot(nr,nc,[5 6]), hold on,
        plot(observ_sig,'color',[0,0,0]+0.45,'linewidth', 1.5),
        text(0,0.8,'(d)','fontsize',20),
    hold off
    legend('Recorded CO$_2$','location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    axis([1 sig_size -0.01 1.2*max(observ_sig(:))]);
    
    subplot(nr,nc,[7 8]), hold on,
        plot(recon_sig,'color',[0, 0.4470, 0.7410],'linewidth', 1.2),
        text(0,1.45,'(e)','fontsize',20),
    hold off
    legend(['Corrected CO$_2$ (Delay: ',num2str(delay_bd+Params.delay),')'],'location','northwest','fontsize',15,'interpreter','latex')
    set(gca,'xtick',[])
    axis([1 sig_size -0.01 1.2*max(recon_sig(:))]);
    
    subplot(nr,nc,[9 10 11 12]), hold on,
        plot(observ_sig,'color',[0,0,0]+0.45,'linewidth', 1.5),
        plot(Abdo+0.8,'color',[0.8500, 0.3250, 0.0980],'linewidth', 1.5),
        plot(recon_sig,'color',[0, 0.4470, 0.7410],'linewidth', 1.2),
        text(0,1.6,'(f)','fontsize',20),
    hold off
    set(gca,'ytick',[])
    set(gca,'xtick',[0 sig_size/2 sig_size],'xticklabel',{'0','200','400(s)'},'fontsize',15)
    
%     subplot(nr,nc,[11 12]), hold on,
%         plot(params.b,'color',[0,0,0]+0.45,'linewidth', 1.5),
%         plot(Abdo+1,'r','linewidth', 1.5),
%         plot(recon_sig,'b','linewidth', 1.2),
%     hold off
%     set(gca,'ytick',[])
    
set(gcf,'Position',[50 100 800 1000])