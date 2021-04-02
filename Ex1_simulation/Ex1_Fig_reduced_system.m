%% Figure 1: Reduced system for respirometry inverse problem with impulse response function

n= 512;
ht = params.ht;
delay = params.delay;
supp = params.supp/5*2;

figure, numRecsDown = 1; numRecsAcross = 3;
i=1; j=1;
ll = 0.08;
subplot('Position',[0.05+(j-1)*1/numRecsAcross 0.09+(numRecsDown-i)*1/numRecsDown .8*1/numRecsAcross .82*1/numRecsDown])
delayadd = 20;
hold on;
    p1=plot([zeros(delayadd,1); ht(1:n/2-delayadd)],'k','linewidth',2);
    p2=plot(ones(100,1),linspace(0,0.3,100),'--k','linewidth',1);
    p3=plot((delay+delayadd)*ones(100,1),linspace(0,ll*1.2,100),'--','linewidth',1,'Color',p2.Color);
    p4=plot((delay+delayadd+supp)*ones(100,1),linspace(0,ll*1.2,100),'--','linewidth',1,'Color',p2.Color);
    p5=plot(linspace(1,(delay+delayadd),100),ll*ones(100,1),'linewidth',1,'Color',p2.Color);
    p6=plot(linspace((delay+delayadd),(delay+delayadd+supp),100),ll*ones(100,1),'linewidth',1,'Color',p2.Color);
    plot([1 3],[ll ll*33/32],'linewidth',1,'Color',p2.Color)
    plot([1 3],[ll ll*31/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd delay+delayadd-2],[ll ll*33/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd delay+delayadd-2],[ll ll*31/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd delay+delayadd+2],[ll ll*33/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd delay+delayadd+2],[ll ll*31/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd+supp delay+delayadd+supp-2],[ll ll*33/32],'linewidth',1,'Color',p2.Color)
    plot([delay+delayadd+supp delay+delayadd+supp-2],[ll ll*31/32],'linewidth',1,'Color',p2.Color)
    text((delay+delayadd)/2-10, ll*1.05,'delay','fontsize',10)
    text(delay+delayadd+(supp)/2-16, ll*1.05,'support','fontsize',10)
hold off
set(gca,'xticklabel',[],'yticklabel',[])
ylabel('$\textbf{h(t)}$','fontsize',13,'interpreter','latex')
xlabel('$\textbf{t}$','fontsize',13,'interpreter','latex')
axis([1 n/3.5 0 ll*1.2])

i=1; j=2;
subplot('Position',[(j-1)*1/numRecsAcross (numRecsDown-i)*1/numRecsDown 1/numRecsAcross 1/numRecsDown])
hold on
    fill([1 10 1],[1 1 10],'k','facealpha',.5)
    l1=plot(11*ones(100,1), linspace(1,10,100),'linewidth',5);
    l2=plot(13*ones(100,1), linspace(1,10,100),'linewidth',5);
    text(5,0.4,'$\textbf{H(h)}$','fontsize',13,'interpreter','latex')
    text(10.85,0.4,'$\textbf{x}$','fontsize',13,'interpreter','latex')
    text(12.85,0.4,'$\textbf{b}$','fontsize',13,'interpreter','latex')
hold off
axis off
axis([0 14 0 11])

i=1; j=3;
subplot('Position',[(j-1)*1/numRecsAcross (numRecsDown-i)*1/numRecsDown 1/numRecsAcross 1/numRecsDown])
hold on
    plot(ones(100,1), linspace(8,10,100),':k','linewidth',1.5)
    plot(ones(100,1), linspace(1,4,100),':k','linewidth',1.5)
    plot(linspace(8,10,100),ones(100,1),':k','linewidth',1.5)
    plot(linspace(1,4,100),ones(100,1),':k','linewidth',1.5)
    plot(linspace(1,10,100),linspace(10,1,100),':k','linewidth',1.5)
    fill([4 8 1 1],[1 1 8 4],'k','facealpha',.5)
    l3=plot(11*ones(100,1), linspace(3,10,100),'linewidth',5);
    plot(11*ones(100,1), linspace(1,3,100),':','linewidth',1.5,'Color',l1.Color);
    l4=plot(13*ones(100,1), linspace(1,8,100),'linewidth',5,'color',l2.Color);
    plot(13*ones(100,1), linspace(8,10,100),':','linewidth',1.5,'color',l2.Color);
    text(5,0.4,'$\widetilde{\textbf{H}}(\widetilde{\textbf{h}})$','fontsize',13,'interpreter','latex')
    text(10.85,0.4,'$\widetilde{\textbf{x}}$','fontsize',13,'interpreter','latex')
    text(12.85,0.4,'$\widetilde{\textbf{b}}$','fontsize',13,'interpreter','latex')
hold off
axis off
axis([0 14 0 11])

set(gcf,'Position',[300 300 600 180])
