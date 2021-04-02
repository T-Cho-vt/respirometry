function [X,INFO] = linsolve_respirometry(H,b,lin_method,lin_method_options)
%
% LINSOLVE_RESPIROMETRY solves linear inverse problem such that
%
%               b = H(h)*x+noise
%
% where - b is output(observation), 
%       - h is a known impulse response,
%       - H(h) is a known forward operator(linear) defined by h,
%       - x is input(original signal),
%       - noise ~ N(0,sigma^-2*I),
%       - sigma^2 is a noise level.
%
% *** Note the length of b is equal to the length of x. ***
%
% X = LINSOLVE_RESPIROMETRY(H,b) solves linear system for X.
% 
% X = LINSOLVE_RESPIROMETRY(H,b,lin_method) solves linear system by a
% chosen method with its default option. The possible methods are:
%               
% X = LINSOLVE_RESPIROMETRY(H,b,lin_method_option) solves lineat system
% by a chosen method with user defined option. 
%
% [X,INFO] = LINSOLVE_RESPIROMETRY(H,b,...) returns output if the chose
% methods provide output other than X.
%
% To use this methods, H and b are required and the rest of them are
% optional. For lin_method and lin_method_options are optional. Unless user 
% defines them, default options are applied. {} indicates default for each 
% option.
%
% Inputs: 
%  H : either (a) a full or sparse matrix
%             (b) a matrix object that performs the matrix*vector operation
%             (c) user-defined function handle
%
%  b : right-hand side vector
%
%  lin_method: solving ill-posed linear inverse problem. (optional)
%
%              [ '\' | {'HyBR'} | 'genHyBR' | 'FLSQR' | 'FISTA']
%
%              - \ is backslash and solving normal equation from the given 
%                linear inverse problem for given either regularization
%                parameter or x_true. Using 2-norm regularization.
%              - HyBR is hybrid iterative projection methods with svd-based
%                regularization method. Use prior I as covariance matrix.
%                Using 2-norm regularization.
%              - genHyBR is extended hybrid methods of HyBR. Use prior
%                covariance matrix Q. Using 2-norm regularization.
%              - FLSQR is hybrid iterative projection methods enforcing
%                1-norm regularization.
%              - FISTA is regularization method to enforce sparsity in the
%                computed X. Using 1-norm regularization.
%
%  lin_method_options: structur with the following fields (optional)
%                      See each method to give options.
%
%       x_true     - true solution; allows us to returns error norms with
%                    respect to x_true 
%       lambda     - parameter for 1-norm (FISTA) or 2-norm (\)
%                    regularization; default is 0.01;
%       L          - weight matrix for X; default = identiy matrix
%       mu         - prior mean of X; default = zero vector
%       Q          - prior covariance matrix; default = identity matrix
%       R          - covariance matrix of noise; default = identity matrix
%       HyBRset    - Options for HyBR. See more details in HyBR.m
%       genHyBRset - Options for genHyBR. See more details in genHyBR.m
%       IRset      - Options for FLSQR. See more details in IRhybrid_flsqr.m
%       FISTA_opts - Options for FISTA. See more details in fista_lasso.m
%
%
% For lin_method,
%      If '\' is chosen (not recommended for too large-scale), 
%         - H must be full matrix,
%         - lin_method_options.L must be given,
%         - either lin_method_options.lambda  or 
%                  lin_method_options.x_true  must be given.
%      If 'HyBR' is chosen,
%         - lin_method_options.mu is optional but recommmend,
%         - lin_method_options.HyBRset must be given. 
%      If 'genHyBR' is chosen,
%         - lin_method_options.mu is optional but recommmend,
%         - lin_method_options.Q must be given,
%         - lin_method_options.R must be given,
%         - lin_method_options.genHyBRset must be given.
%      If 'FLSQR' is chosen,
%         - lin_method_options.mu is optional but recommmend,
%         - lin_method_options.IRset must be given. 
%      If 'FISTA' is chosen,
%         - H must be full matrix,
%         - lin_method_options.lambda must be given,
%         - lin_method_options.FISTA_opts must be given.
%
% Outputs:
%  X : computed solutions
%  INFO : more information of output
%
% Taewon Cho, Virgina Tech
% October, 2020.

n = length(b);

if nargin < 4
    
    lin_method_options.x_true = [];
    lin_method_options.lambda = [];
    lin_method_options.mu = [];
    lin_method_options.L = [];
    lin_method_options.Q = [];
    lin_method_options.R = [];
    lin_method_options.HyBRset = [];
    lin_method_options.genHyBRset = [];
    lin_method_options.IRset = [];
    
    if nargin < 3
        lin_method = 'HyBR';
    end
end

if isempty(lin_method)
   lin_method = 'HyBR';
end

if ~myIsField(lin_method_options,'x_true')
    lin_method_options.x_true = [];
end
if ~myIsField(lin_method_options,'lambda')
    lin_method_options.lambda = [];
end
if ~myIsField(lin_method_options,'mu')
    lin_method_options.mu = [];
end
if ~myIsField(lin_method_options,'L')
    lin_method_options.L = [];
end
if ~myIsField(lin_method_options,'Q')
    lin_method_options.Q = [];
end
if ~myIsField(lin_method_options,'R')
    lin_method_options.R = [];
end
if ~myIsField(lin_method_options,'HyBRset')
    lin_method_options.HyBRset = [];
end
if ~myIsField(lin_method_options,'genHyBRset')
    lin_method_options.genHyBRset = [];
end
if ~myIsField(lin_method_options,'IRset')
    lin_method_options.IRset = [];
end


switch lin_method
    
    case '\'
        % Check weight matrix for X
        if isempty(lin_method_options.L)
            lin_method_options.L = speye(n);
        end
        % Check if x_true exists
        if isempty(lin_method_options.x_true)
            lin_method_options.x_true = zeros(n,1);
            % Check if lambda is given
            if isempty(lin_method_options.lambda)
                lin_method_options.lambda = 0.01;
            end
            % H must be a full matrix.
            X = (H'*H+lin_method_options.lambda^2*lin_method_options.L'*lin_method_options.L)\(H'*b);
            if nargout >1
                INFO.lambda = lin_method_options.lambda;
            end
        else
            % H must be a full matrix.
            gtik = @(lambda) norm(((H'*H+lambda^2*L'*L)\(H'*b)) - lin_method_options.x_true);
            tik_lambda = fminbnd(gtik,0,1);
            f_Tik_Q = (H'*H+tik_lambda^2*L'*L)\(H'*b);
            X = f_Tik_Q;
            if nargout >1
                INFO.lambda = tik_lambda;
            end
        end
        
    case 'HyBR'
        % Check the prior mean
        if isempty(lin_method_options.mu)
            lin_method_options.mu = zeros(n,1);
        end
        % Check the options for HyBR
        if isempty(lin_method_options.HyBRset)
            input_hybr = HyBRset('HyBR');
            % Check if x_true exists
            if isempty(lin_method_options.x_true)
                lin_method_options.HyBRset = HyBRset(input_hybr,'Iter',100);
            else
                lin_method_options.HyBRset = HyBRset(input_hybr,'Iter',100,'x_true',lin_method_options.x_true-lin_method_options.mu);
            end
        end
        d = b - H*lin_method_options.mu;
        [Y,output] = HyBR(H,d,[],lin_method_options.HyBRset);
        X = Y + lin_method_options.mu;
        if nargout >1
            INFO = output;
        end
        
    case 'genHyBR'
        % Check the prior mean
        if isempty(lin_method_options.mu)
            lin_method_options.mu = zeros(n,1);
        end
        % Check the options for genHyBR
        if isempty(lin_method_options.genHyBRset)
            input_genhybr = HyBRset('genHyBR');
            % Check if x_true exists
            if isempty(lin_method_options.x_true)
                lin_method_options.genHyBRset = HyBRset(input_genhybr,'Iter',100);
            else
                lin_method_options.genHyBRset = HyBRset(input_genhybr,'Iter',100,lin_method_options.x_true-lin_method_options.mu);
            end
        end
        % Check the prior covariance matrix
        if isempty(lin_method_options.Q)
            lin_method_options.Q = speye(n);
        end
        % Check the covariance matrix of noise
        if isempty(lin_method_options.R)
            lin_method_options.R = speye(n);
        end
        d = b - H*lin_method_options.mu;
        [Y,output] = genHyBR(H,d,lin_method_options.Q,lin_method_options.R,lin_method_options.genHyBRset);
        X = Y + lin_method_options.mu;
        if nargout >1
            INFO = output;
        end
        
    case 'FLSQR'
        % Check the prior mean
        if isempty(lin_method_options.mu)
            lin_method_options.mu = zeros(n,1);
        end
        % Check the options for FLSQR
        if isempty(lin_method_options.IRset)
            input_flsqr = IRset('IRhybrid_flsqr');
            % Check if x_true exists
            if isempty(lin_method_options.x_true)
                lin_method_options.IRset = IRset(input_flsqr,'MaxIter',100,'hybridvariant','R');
            else
                lin_method_options.IRset = IRset(input_flsqr,'MaxIter',100,'hybridvariant','R','x_true',lin_method_options.x_true);
            end
        end
        d = b - H*lin_method_options.mu;
        [Y,output] = IRhybrid_flsqr(H,d,[],lin_method_options.IRset);
        X = Y + lin_method_options.mu;
        if nargout >1
            INFO = output;
        end
        
    case 'FISTA'
        % Check the prior mean
        if isempty(lin_method_options.mu)
            lin_method_options.mu = zeros(n,1);
        end
        
        d = b - H*lin_method_options.mu;
        
        % Check if lambda is given
        if isempty(lin_method_options.lambda)
            [Y,INFO_fista]=IRfista(H,d);
        else
            options=IRset('RegParam',lin_method_options.lambda);
            [Y,INFO_fista]=IRfista(H,d,[],options);
        end
        
        X = Y + lin_method_options.mu;        
        
        if nargout >1
            INFO = INFO_fista;
        end
            
end

