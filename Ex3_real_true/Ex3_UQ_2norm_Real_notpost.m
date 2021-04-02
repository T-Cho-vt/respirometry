function UQ = Ex3_UQ_2norm_Real(output,nlevel)
% UQ_2-norm
% Tikhonov
UQ.lambda = 0.1;
UQ.nLevel = nlevel;
UQ.A_psue = inv(output.Hfull'*output.Hfull+UQ.lambda^2*output.L'*output.L)*output.Hfull';
UQ.x_tik = UQ.A_psue*output.b;
UQ.Var = (UQ.nLevel^2)*(UQ.A_psue*UQ.A_psue');
clear UQ.A_psue;
UQ.x_std = sqrt(abs(diag(UQ.Var)));
clear UQ.Var;

min_fig = min( [min(UQ.x_tik), min(output.b)] );
max_fig = max( [max(UQ.x_tik), max(output.b)] );

% Figure
UQ.color1  = [0 0 0];
UQ.color2 = [0, 0.4470, 0.7410];
figure, hold on
    A1 = UQ.x_tik+1.96*UQ.x_std; 
    A2 = UQ.x_tik-1.96*UQ.x_std;   
    h1 = fill([1:1:output.n,output.n:-1:1],[A1', fliplr(A2')],UQ.color1);
    set(h1,'facealpha',0.2,'EdgeColor',UQ.color1)
    h2 = plot(output.b,'linewidth',2,'color',UQ.color2);
    h3 = plot(1:1:output.n,UQ.x_tik,'r','linewidth',1);
    h4=fill([1 1 output.delay output.delay], [min_fig max_fig max_fig min_fig] , 'r',...
        [output.n-output.delay+1 output.n-output.delay+1 output.n output.n],...
        [min_fig max_fig max_fig min_fig] , 'r');
    set(h4,'FaceAlpha',0.1);
    set(h4,'Edgecolor','none');
    h5=fill([output.delay+1 output.delay+1 output.delay+output.supp output.delay+output.supp], ...
        [min_fig max_fig max_fig min_fig] , 'b',...
        [output.n-output.supp-output.delay+1 output.n-output.supp-output.delay+1 output.n-output.delay+1 output.n-output.delay+1], ...
        [min_fig max_fig max_fig min_fig] , 'b');
    set(h5,'FaceAlpha',0.1);
    set(h5,'Edgecolor','none');
hold off
% axis([1 n 0 1])
h6=legend([h2 h3 h1],{'Observation','Tikhonov','95% c.b.'},'Orientation','horizontal');
title([num2str(100*UQ.nLevel),'% noise'])
set(gcf,'Position',[300 300 1200 400])