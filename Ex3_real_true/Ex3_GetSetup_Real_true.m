function output = GetSetup_Real_true(obsize)

SamplingRate = 100;
% experimental data
P = textread('y.txt');
P=P(1:obsize,:); % selecting a part of the data; obsize 20000 or 40000
P = P(1:100/SamplingRate:end,:);    % sampling rate has been 100 Hz. We are making it 10Hz here
output.b = P(:);

n = length(output.b);
output.n = n;

% Impulse response
h=textread('h.txt');
h=h(1:100/SamplingRate:end)';    % making the sampling rate 10Hz
ht=h/sum(h)*SamplingRate;
output.ht = zeros(n, 1);
output.ht(1:length(ht)) = ht/SamplingRate;
output.ht = output.ht(:);

% true signal
x_exact = textread('u_exact.txt');
output.x_exact = x_exact(1:obsize);

% -------- Define H matrix --------
output.bndry = 'reflexive'; % This can be {'zero'|'periodic'|'reflexive'}
% Hfull = matrixH(output.ht,output.bndry);
% output.Hfull = Hfull;
H = matrixHfun(output.ht,output.bndry,1);
output.H = H;

% -------- Get delay and support of true h
spt = output.ht>0.0001; 
% spt = output.ht>0.001; 
% spt = output.ht>0.005; 
fndnz = find(spt==1); 
delay = fndnz(1)-1; % delay
output.delay = delay;
supp = sum(spt); % number of entries in support
output.supp = supp;

% -------- Weight matrix for x in penalty norm
% first derivative operator
Q1 = Q1_Matrix(n);
% second derivative operator
Q2 = Q2_Matrix(n);      
% Identity
I = speye(n);
% choose either {'Q1'|'Q2'|'I'}
output.L = Q1;

end