%% Tikhonov regularization for x 
Params.solve = 'HyBR';
Params.MaxIter = 0;
Params.finallin = 'off';
h0 = ones(params.supp,1)/params.supp; % initial guess
[mot_x,mot_h,mot_res,mot_flag] = Ex1_blind_sim(h0,Params);

% Plot x, h, and b
figure, 
numRecsDown = 2;
numRecsAcross = 3;

i=1; j=1;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 1.1*(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
plot(params.x_exact,'linewidth',1.5)
axis([1 length(params.x_exact) -0.1 1.1])
title('x_{true}')

i=1; j=2;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 1.1*(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
plot(params.ht,'linewidth',1.5)
title('h_{true}')
axis([1 n 0 1.5*max(params.ht)])

i=1; j=3;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 1.1*(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
plot(params.b,'linewidth',1.5)
title('b_{true}')
axis([1 n -0.2 1])

colNum = 1;
h_aug = [zeros(params.delay,1);mot_h; zeros(n-params.delay-params.supp,1)];
Hcheck = matrixH(h_aug,Params.bndry);
diff_b = norm( Hcheck*mot_x- params.b)/norm(params.b);

i=2; j=3;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 0.05+(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
hold on
plot(params.b,'linewidth',1.5)
p1 = plot(Hcheck*mot_x,'--','linewidth',1.5);
hold off
legend('b_{true}','b_1','location','northwest')
title(['\Delta b_1 = ',sprintf('%.1d',diff_b)])
axis([1 n -0.2 1])

i=2; j=1;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 0.05+(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
plot(mot_x,'color',p1.Color,'linewidth',1.5)
title('x_1')
axis([1 n 1.1*min(mot_x), 1.1*max(mot_x)])

i=2; j=2;
subplot('Position',[0.04+(j-1)*1/numRecsAcross 0.05+(numRecsDown-i)*1/numRecsDown 0.8*1/numRecsAcross 0.8*1/numRecsDown])
hold on
plot(h_aug,'color',p1.Color,'linewidth',1.5), hold off
axis([1 length(h_aug) 0 0.25])

title('h_1')
axis([1 n 0 1.5*max(params.ht)])
set(gcf,'Position',[300 300 600 400])