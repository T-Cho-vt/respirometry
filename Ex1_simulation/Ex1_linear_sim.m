function output_lin = Ex1_linear_sim(n,output)
% Test linear system solvers

%% Tikhonov - Q
gtik = @(alpha) norm( ( (output.Hfull'*output.Hfull+alpha^2*output.Q'*output.Q)\(output.Hfull'*output.b))- output.x_exact);
tik_alpha = fminbnd(gtik,0,1);
f_Tik_Q = (output.Hfull'*output.Hfull+tik_alpha^2*output.Q'*output.Q)\(output.Hfull'*output.b);
output_lin.f_Tik_Q = f_Tik_Q;
output_lin.alpha_Tik_Q = tik_alpha;

figure,
subplot(4,1,1), 
pb1=plot(output.x_exact,'linewidth', 1.5); hold on, pb2=plot(f_Tik_Q,'linewidth', 1.5); hold off, 
title(['Tikhonov-Q-opt, Err = ',num2str(norm(output.x_exact-f_Tik_Q)/norm(output.x_exact))])
axis([1 n -0.2 1.3])


%% HyBR 
mu=zeros(n,1);
d = output.b-output.H*mu;
input_hybr = HyBRset('HyBR');
input_hybr = HyBRset(input_hybr,'Iter',100,'x_true',output.x_exact);
[y,output_HyBR] = HyBR(output.H,d,[],input_hybr);
f_HyBR = y + mu;
output_lin.f_HyBR = f_HyBR;
output_lin.output_HyBR = output_HyBR;

% figure,
subplot(4,1,2), 
pb1=plot(output.x_exact,'linewidth', 1.5); hold on, pb2=plot(f_HyBR,'linewidth', 1.5); hold off, 
title(['HyBR-I, Err = ',num2str(norm(output.x_exact-f_HyBR)/norm(output.x_exact))])
axis([1 n -0.2 1.3])
% set(gcf,'Position',[300 300 600 120])

%% FLSQR
FLSQR_opt.x_true = output.x_exact;
FLSQR_opt = IRset(FLSQR_opt, 'hybridvariant', 'R');% ,'SparsityTrans','dwt');
[f_FLSQR,output_FLSQR] = IRhybrid_flsqr(output.H,output.b,[],FLSQR_opt);
output_lin.f_FLSQR = f_FLSQR;
output_lin.output_FLSQR = output_FLSQR;

% figure, 
subplot(4,1,3),
plot(output.x_exact,'linewidth', 1.5), hold on,
plot(f_FLSQR,'linewidth', 1.5), hold off,
title(['FLSQR-R, Err = ',num2str(norm(output.x_exact-f_FLSQR)/norm(output.x_exact))])
axis([1 n -0.2 1.3])

%% FISTA
opts.lambda = 0.0005;
f_fista_lasso = fista_lasso(output.b,output.Hfull,[], opts);
output_lin.f_fista_lasso = f_fista_lasso;
output_lin.lambda_fista = opts.lambda;

% figure,
subplot(4,1,4),
plot(output.x_exact,'linewidth', 1.5), hold on, plot(f_fista_lasso,'linewidth', 1.5), hold off, 
title(['FISTA, Err = ',num2str(norm(output.x_exact-f_fista_lasso)/norm(output.x_exact))])
axis([1 n -0.2 1.3])