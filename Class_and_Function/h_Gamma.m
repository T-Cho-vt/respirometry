function [H, h, h_dot] = h_Gamma(y, T, f, tol, L)

%  h    : generated inpulse response,   h= h(alpha,beta) 
%  h_dot: derivative of the generated impulse response wrt alpha and beta
%  T    : duration of the impulse response (like 50-100 seconds)
%  f    : sampling frequency (unusally use 10 Hz)
%  tol  : tolerance to eliminate the first zeros (or close to zero) of h
%  L    : length of h. This is not for gernerating the true observed signal.
%         Before running nonlinear solver, this function should be set up 
%         as blur function with L

    Alpha = y(1); Beta = y(2);

    t = (1/f:1/f:T)';
    scl = (Beta^(Alpha+1))/gamma(Alpha+1);
    h = scl * t.^Alpha .* exp(-Beta*t);
    
    % Cut zeros (less than tol) at the begining of h
    if nargin < 5
        i = 1;
        while h(i) < tol
            i = i+1 ;
        end
        h = h(i:end);  
    elseif nargin == 5 % Match length of h between h_0 and h_k
        i = length(h)-L+1;
        h = h(i:end);
    end
    
    % Construct H explictly (not efficient)
    n = length(h);
    H = toeplitz(h, [h(1);zeros(n-1,1)])/f; 

    if nargout > 2
        h_dot = [(log(Beta) + log(t(i:end)) - psi(Alpha+1)) .* h, ...
                 ((Alpha+1)/Beta - t(i:end)) .* h];
    end
  
end
