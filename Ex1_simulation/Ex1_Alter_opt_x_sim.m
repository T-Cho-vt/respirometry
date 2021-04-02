function [X_red, x] = Ex1_Alter_opt_x_sim(h0, param)
%
% Input:
%       h0 - impulse function at given iteration on support
%    param - has the following fields
%                     b : observed signal
%                nLevel : noise level
%                     T : true signal
%                  beta : regularization parameter for support of x for '\'
%                 delay : delay of h
%                  supp : length of support
%                 solve : (linear) solver method to compute signal, x
%               MaxIter : max iterations for outer iteration
%                 solve : (linear) solver method to compute signal, x
%
% Output:
%    X_red - Reduced cost functional
%        x - singal for given h0
%
% T. Cho 5/21/2020

b = param.b;
n = size(b,1);
delay = param.delay;
supp = param.supp;
% useless observation for delay
b = b(param.delay+1:end);

% Generate forward operator H for given h0.
h0_aug = [h0; zeros(n-delay-supp,1)];

H = matrixH(h0_aug,param.bndry);

% Solve for x
switch param.solve
        
    case 'FISTA'
        opts.lambda = param.beta; 
        x = fista_lasso(b,H,[],opts);
    
    case 'FLSQR'
        opt.x_true = param.T(1:n-delay);
        opt = IRset(opt, 'hybridvariant', 'R');
        [x,x_info] = IRhybrid_flsqr(H,b,[],opt);
        
    case 'HyBR'
        input_hybr = HyBRset('HyBR');
        input_hybr = HyBRset(input_hybr,'Iter',100);
        [x,x_output] = HyBR(H,b,[],input_hybr);
        
end

% Compute reduced function value
Hx = H*x(:);
X_red = 0.5*(Hx - b)'*(Hx - b);
if nargout >1
    % Get full x
    x = [x; zeros(delay,1)];
end