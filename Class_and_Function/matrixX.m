function X = matrixX(x,cond)

x=x(:);
n = length(x);

switch cond
    case 'zero'
        X = toeplitz(x, [x(1);zeros(n-1,1)]); 
    case 'periodic'
        X = toeplitz(x, [x(1);flipud(x(2:end))]); 
    case 'reflexive'
        X = toeplitz(x, [x(1);x(1:end-1)]); 
end