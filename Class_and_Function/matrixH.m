function H = matrixH(h,cond)

h = h(:);
n = length(h);

switch cond
    case 'zero'
        H = toeplitz(h, [h(1);zeros(n-1,1)]); 
    case 'periodic'
        H = toeplitz(h, [h(1); flipud(h(2:end))]);
    case 'reflexive'
        H = toeplitz(h, [h(1);zeros(n-1,1)]) + fliplr(toeplitz(zeros(n,1), [0; flipud(h(2:end))]));
end
