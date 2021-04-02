function [X_red, x] = Ex2_Alter_opt_x_Real(h0, param)
%
% Input:
%       h0 - impulse function at given iteration on support
%    param - has the following fields
%                     b : observed signal
%                nLevel : noise level
%                     T : true signal
%                  beta : regularization parameter for x
%                 delay : delay of h
%                  supp : length of support
%                 solve : (linear) solver method to compute signal, x
%               MaxIter : max iterations for outer iteration
%
% Output:
%    X_red - Reduced cost functional
%        x - singal for given h0
%
% T. Cho 5/21/2020

b = param.b;
n = size(b,1);
delay = param.delay;

delay = 0;

supp = param.supp;
% useless observation for delay
b = b(delay+1:end);

% Generate forward operator H for given h0.
h0_aug = [h0; zeros(n-delay-supp,1)];
h0_aug = [zeros(param.delay,1); h0; zeros(n-param.delay-supp,1)];

% H = matrixH(h0_aug,param.bndry);
H = matrixHfun(h0_aug,param.bndry,1);

% Solve for x
switch param.solve
        
    case 'FISTA'
        opts_fista.lambda = param.beta;
        opts_fista.mu = zeros(n,1);
        x = linsolve_respirometry(H,b,'FISTA',opts_fista);
    
    case 'FLSQR'
         [x,x_output] = linsolve_respirometry(H,b,'FLSQR');
        
    case 'HyBR'
         opts_hybr.mu = sum(b)/n *ones(n,1); 
         [x,x_output] = linsolve_respirometry(H,b,'HyBR',opts_hybr);
end

% Compute reduced function value in Variable Projection
Hx = H*x(:);
X_red = 0.5*(Hx - b)'*(Hx - b);
if nargout >1
    % Get full x
    x = [x; zeros(delay,1)];
end