function [x, h0, Hs, reserr, flag] = Ex1_blind_sim(h0,param)
%
% Solve separable nonlinear problem with variable projection
%
% Input:
%       h0 - initial guess for supported h
%    param - has the following fields
%                     b : observed signal
%                nLevel : noise level
%                     T : true signal
%                    ht : true h (or experimental h)
%                 alpha : regularization parameter for support of h
%                  beta : regularization parameter for support of x
%                 delay : delay of h
%                  supp : length of support of h
%                 solve : (linear) solver method to compute signal, x
%           nonlinsolve : nonlinear solver for h
%              RegMat_h : weight matrix for support of h 
%              RegMat_x : weight matrix for support of x 
%                  hnrm : h norm change stopping
%                  rnrm : residual norm error change stopping
%               MaxIter : max iterations for outer iteration
%              finallin : Do one more linear optimization from the final h
%                         {'on'|'off'}
%              plotshow : Show x and h while running codes
%
% Output:
%        x - final computed signal(image)
%       h0 - final computed h
%       Hs - store all h in matrix H0
%   reserr - residual errors
%     flag - reason for stopping iteration
%             1: change of h_i < param.hnrm
%             2: change of residual < param.rnrm
%             3: reach the maximum iteration
%
% T.Cho 05/21/2020


b = param.b;
n = size(b,1); % size of signal(f) and b
supp = param.supp; 

% True h
ht = param.ht;

% Stack gamma and signals, residual error
Hs(:,1) = h0(:);
Signals = [];
reserr = [];

% Iteration count
iter = 0; 
resold = 0;
if param.plotshow
    figure;
end

while 1
    % Update count
    iter = iter+1;
    param.iter = iter; % to check iteration in Ex1_Alter_opt_x_sim
    % Evaluate signal for given h0
    [~,x] = Ex1_Alter_opt_x_sim(h0, param);
    param.f_b4 = x; % store x for next iteration
    
    % Store x at each iteration
    Signals(:,iter) = x(:);
    
    if param.plotshow
        % Plot the current x,h
        % plot x
        str = ['itr = ',int2str(iter)];
        subplot(1,2,1), 
            plot(param.T(1:end-param.delay)), 
            hold on, 
                plot(x(1:end-param.delay)), 
            hold off;
            xlim([1 length(param.b(1:end-param.delay))])
            title(['Signal, ',str]);
    end
    
    % Choose nonlinear solver
    switch param.nonlinsolve
            
        case 'lsqlin'
            % Use lsqlin
            Toe = matrixH(x(1:end-param.delay),param.bndry);
            Toe_supp = Toe(:,1:supp);
            Toeb = param.b(param.delay+1:end);
            
            % Add reg
            Toe_supp = [Toe_supp; param.alpha*eye(supp)];
            Toeb = [Toeb; zeros(supp,1)];
            
            hest = lsqlin(Toe_supp, Toeb, [], [],...
                [ones(1,supp); [1 zeros(1,supp-1)];[ zeros(1,supp-1) 1]],...
                [1;0;0],...
                zeros(supp,1),inf*ones(supp,1));
    end

    % Update fields
    hOld = h0;
    h0 = hest;
    Hs(:,iter+1) = h0(:);
    
    if param.plotshow
        % plot h
        subplot(1,2,2)
            plot(ht(param.delay+1:param.delay+param.supp))
            hold on,
                plot(h0)
            hold off; 
            title(['h, ',str]);
            xlim([1 param.supp])
            pause(0.01);
    end

    % Check for termination: Note two criteria.
    %           (1) norm(hOld(:)-h0(:),'inf') < 1e-6
    %           (2) abs(resold-resnew) < 1e-6
    %           (3) iter > param.MaxIter
    
    if norm(hOld(:)-h0(:)) < param.hnrm
        flag=1; 
        switch param.finallin
            case 'on'
                [x,h0] = searching_reg_parameter(h0,param);
                resh = [zeros(param.delay,1); h0; zeros(n-param.delay-supp,1);];
                resH = matrixHfun(resh,1);
                resnew = norm(resH*x-b);
                reserr = [reserr; resnew];
                return; 
            case 'off'
                return; 
        end
    end  
    
    resh = [zeros(param.delay,1); h0; zeros(n-param.delay-supp,1);];
    resH = matrixHfun(resh,param.bndry,1);
    resnew = norm(resH*x-b);
    reserr = [reserr; resnew];
	
    if abs(resold-resnew) < param.rnrm 
        flag=2; 
        switch param.finallin
            case 'on'
                [x,h0] = searching_reg_parameter(h0,param);
                resh = [zeros(param.delay,1); h0; zeros(n-param.delay-supp,1);];
                resH = matrixHfun(resh,1);
                resnew = norm(resH*x-b);
                reserr = [reserr; resnew];
                return; 
            case 'off'
                return; 
        end
    end
    resold = resnew;
    
    if iter > param.MaxIter 
        flag=3; 
        switch param.finallin
            case 'on'
                [x,h0] = searching_reg_parameter(h0,param);
                resh = [zeros(param.delay,1); h0; zeros(n-param.delay-supp,1);];
                resH = matrixHfun(resh,1);
                resnew = norm(resH*x-b);
                reserr = [reserr; resnew];
                return; 
            case 'off'
                return; 
        end
    end
end