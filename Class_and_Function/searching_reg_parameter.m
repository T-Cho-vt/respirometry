function [x,h0] = searching_reg_parameter(h0,param)


    H = matrixH([zeros(param.delay,1); ...
                 h0; ...
                 zeros(length(param.ht)-param.delay-param.supp,1)]);

    fun = @(alpha) BestRegFISTA(alpha,H,param.b,param.T);
    best_fista_alpha = fminbnd(fun,0,0.1);
    opts.lambda = best_fista_alpha;
    x = fista_lasso(param.b,H,[],opts);