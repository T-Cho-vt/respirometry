function output = GetSetup_Real(obsize)

SamplingRate = 10;
% experimental data
P = textread('Data_10Hz.txt');
P=P(1:obsize,:); 
P = P(1:10/SamplingRate:end,:);    % sampling rate has been 10 Hz. We are making it 10Hz here
output.b = P(:,2)-0.25; % shifted observation close to zero but no negative

n = length(output.b);
output.n = n;

% Impulse response
h=textread('ImpulseResponse_100ms.txt');
h=h(1:100/SamplingRate:end)';    % making the sampling rate 10Hz
ht=h/sum(h)*SamplingRate;
output.ht = zeros(n, 1);
output.ht(1:length(ht)) = ht/SamplingRate;
output.ht = output.ht(:);

% -------- Define H matrix --------
output.bndry = 'reflexive'; % This can be {'zero'|'periodic'|'reflexive'}
H = matrixHfun(output.ht,output.bndry,1);
output.H = H;

% -------- Get delay and support of true h
spt = output.ht>0.0001; 
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