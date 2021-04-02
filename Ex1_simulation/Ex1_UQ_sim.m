function output_UQ = Ex1_UQ_sim(n,output)

%% Get UQ w/ 2-norm

gtik = @(alpha) norm( ( (output.Hfull'*output.Hfull+alpha^2*output.Q'*output.Q)\(output.Hfull'*output.b))- output.x_exact);
alpha = fminbnd(gtik,0,1);
output_UQ.UQ_alpha = alpha;
A_psue = inv(output.Hfull'*output.Hfull+alpha^2*output.Q'*output.Q)*output.Hfull';
x_tik = A_psue*output.b;
output_UQ.UQ_x_tik = x_tik;
Var = (output.nLevel^2)*(A_psue*A_psue');
x_std = sqrt(abs(diag(Var)));
output_UQ.UQ_x_std = x_std;

figure, 
hold on
    A1 = x_tik+1.96*x_std; A2 = x_tik-1.96*x_std;
    color  = [0 0 0];
    color2 = [0, 0.4470, 0.7410];
    h = fill([1:1:n,n:-1:1],[A1', fliplr(A2')],color);
    set(h,'facealpha',0.2,'EdgeColor',color)
    h2 = plot(output.x_exact,'linewidth',2,'color',color2);
    h3 = plot(1:1:n,x_tik,'r','linewidth',1);
hold off
axis([1 n -0.4 1.3])
h4=legend([h2 h3 h],{'true image','Tikhonov - Q','95% c.b.'},'location','northoutside','Orientation','horizontal');
set(gcf,'Position',[300 300 600 200])

%% Get UQ w/ 1-norm approx
N=128;
nLevel = 0.005;
output_uq1 = Ex1_GetSetup_sim(N,nLevel);

t = linspace(1,N,N)';
sigma = nLevel * norm(output_uq1.Hfull*output_uq1.x_exact) / sqrt(N);

%  Haar matrix for regularization function and regularization parameter.
D                = HaarWaveletTransformMatrix(N);
D = eye(N);
Nresid           = N + length(D(:,1));
% Now compute the MAP estimator using the nonlinear transformation approach
p.A              = output_uq1.Hfull;
p.lam            = 1/sigma^2;
p.del            = 2;
p.b              = output_uq1.b;
p.D              = D;
p.Q              = speye(Nresid);     % this is non-identity when we do RTO.
p.e              = zeros(Nresid,1);   % this is iid Gaussian when we do RTO.
u0               = ones(N,1); % initial guess
[uMAP,rMAP,JMAP] = LevMar(u0,@(u)l1_fun(u,p),0.001,1e-8,100);
% mse              = norm(rMAP)^2/(M-N);
xMAP             = D'*g_fn(uMAP,p.del);

% Plot the measurements and model fit.
% figure, plot(t,output_uq1.x_exact,t,xMAP)

% Now, use RTO-MH to sample from the posterior distribution defined by
%              
%            p(x|b) \propto exp(-lam/2*||A(x)-b||^2-del/2*||D*x||^2).
%
nsamp                = 1000;
[Q,~]                = qr(JMAP,0);
p.Q                  = Q;
p.Nrand              = Nresid;
[zchain,accept_rate] = RTO_MH(uMAP,@(u,p)l1_fun(u,p),p,nsamp);
xchain               = D'*g_fn(zchain,p.del);

%% Visualize the MCMC chain

% Plot the sample mean and 95% credibility intervals.
xlims          = plims(xchain',[0.025,0.5,0.975])';
relative_error = norm(output_uq1.x_exact-xlims(:,2))/norm(output_uq1.x_exact);
figure,
hold on
plot(t,output_uq1.x_exact,'linewidth',2,'color',color2)
plot(t,xlims(:,2),'r','linewidth',1)
color  = [0 0 0];
h = fill([1:1:N,N:-1:1],[xlims(:,1)', fliplr(xlims(:,3)')],color);
set(h,'facealpha',0.2,'EdgeColor',color)
% axis([1 N 1.1*min(xlims(:,1)) 1.1*max(xlims(:,3))])
axis([1 N -.38 1.304])
set(gcf,'Position',[300 300 600 200])
legend('true image','MCMC sample median','95% credibility bounds','location','northoutside','Orientation','horizontal')
