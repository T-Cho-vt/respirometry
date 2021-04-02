function Params = Ex1_Param_Setup_sim(n,output)

Params.b = output.b; % Observed signal
Params.nLevel = []; % There is no noise level
Params.T = output.x_exact; % There is no true signal
Params.ht = output.ht; % experimental values of h
Params.alpha = 1;%*length(ht); % regularization parameter for h
Params.beta = 0.002; % regularization parameter for x
Params.delay = output.delay; % delay of h
Params.supp = output.supp; % length of support
Params.bndry = 'zero'; % boundary condition of x
Params.plotshow = 0;

% --- Solver for x ---
Params.solve = 'FISTA'; 
% Params.solve = 'FLSQR'; 
% Params.solve = 'HyBR';

% --- Solver for h ---
Params.nonlinsolve = 'lsqlin';

% Use a Regularization (weight) matrix for h
Params.RegMat_h = ones(1,output.supp);
% Use a regularization (weight) matrix for x
Params.RegMat_x = ones(1,n);

% Stopping criteria
Params.hnrm = 1e-6; % Relative error criteria for h
Params.rnrm = 1e-6; % Residual error
Params.MaxIter = 500; % max iteration of nonlinear solver

% On/Off linsolve for the final h
Params.finallin = 'off';

% weight matrix for x
Params.L1 = output.Q;

% weight matrix for h
 % Select Q(Matern) as Kronecker product
 xmin = 0;           %Coordinates of left corner
 xmax = 1;             %Coordinates of right corner
 nvec = output.supp;
 theta = 1;      %For now set them as isotropic
 % Matern kernel
 nu1 = 1; ell1 = .5;
 k1 = @(r) matern(r,nu1,ell1);
 Qr1 = createrow(xmin,xmax,nvec,k1,theta);% figure, subplot(1,2,1), plot(Qr1)
 Q1fun = @(x)toeplitzproduct(x, Qr1, nvec);
 Q1 = funMat(Q1fun,Q1fun, nvec, nvec);
 invQ1 = inv(Q1*eye(output.supp));% subplot(1,2,2), imagesc(Q1*eye(output.supp))
 Params.L2 = chol(invQ1);