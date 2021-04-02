function [X_red, x] = Ex3_Alter_opt_x_Real(h0, param)
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
%              RegMat_x : weight matrix for support of x for '\'
%
% Output:
%    F_red - Reduced cost functional
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
    case '\'
        % param.beta is regularizationn parameter for x. 
        % param.RegMat_x is weight matrix for x
        x = [H; param.beta*diag(param.RegMat_x(delay+1:end))] \ [b; zeros(n-delay,1)];
        % Unconstraint optimization for x causes worst scenario for blind
        % deconvolution problem. x from this method has nonzero tails while
        % the true x has zero tails. For this reason, constraint is
        % strongly recommended -> Use 'lsqlin'
        
    case 'lsqlin'
        % What would be the upper bound for x? Can I give 1 as maximum?
        % When upper bound is inf, h becomes delta function. That's reason
        % why I put 1 as upper bound for x.
        lb = zeros(size(H,2),1);
        ub = 1*ones(size(H,2),1);
        x = lsqlin(H,b, [],[], [],[], lb, ub);
        
    case 'SpaRSA'
        x = SpaRSA(H,b,param.beta);
        x = x';
        
    case 'FISTA'
%         opts.lambda = param.beta; 
%         if param.iter > 1
%             norm2 = norm(param.f_b4);
%             opts.lambda = opts.lambda / norm2;
%         end
%         x = fista_lasso(b,H,[],opts);
        
        opts_fista.lambda = param.beta;
        opts_fista.mu = sum(b)/n *ones(n,1); 
        opts_fista.mu = zeros(n,1);
        x = linsolve_respirometry(H,b,'FISTA',opts_fista);
    
    case 'FLSQR'
%         mu = sum(b(:))/length(b)*ones(n-delay,1);
%         d = b-H*mu;
%         opt.x_true = param.T(1:n-delay);
%         opt = IRset(opt, 'hybridvariant', 'R');
% %         [x,x_info] = IRhybrid_flsqr(H,b,[],opt);
%         [y,y_info] = IRhybrid_flsqr(H,d,[],opt);
%         x = y + mu;
         [x,x_output] = linsolve_respirometry(H,b,'FLSQR');
        
    case 'HyBR'
%         mu = sum(b(:))/length(b)*ones(n-delay,1);
%         d = b-H*mu;
%         input_hybr = HyBRset('HyBR');
%         input_hybr = HyBRset(input_hybr,'Iter',100);
% %         [x,x_output] = HyBR(H,b,[],input_hybr);
%         [y,y_output] = HyBR(H,d,[],input_hybr);
%         x = y + mu;
         opts_hybr.mu = sum(b)/n *ones(n,1); 
         [x,x_output] = linsolve_respirometry(H,b,'HyBR',opts_hybr);
    case 'genHyBR'
        mu = 0.3*ones(n-delay,1);
        d = b-H*mu;
        R = eye(length(b));
        input_genhybr = HyBRset('genHyBR');
        input_genhybr = HyBRset(input_genhybr,'Iter',100);
%         [x,x_output] = genHyBR(H,b,param.Q(1:n-param.delay,1:n-param.delay),R,input_genhybr);
        [y,y_output] = genHyBR(H,d,param.Q(1:n-delay,1:n-delay),R,input_genhybr);
        x = y + mu;
end

% Compute reduced function value in Variable Projection
Hx = H*x(:);
X_red = 0.5*(Hx - b)'*(Hx - b);
if nargout >1
    % Get full x
    x = [x; zeros(delay,1)];
end