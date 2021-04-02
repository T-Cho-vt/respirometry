function output = GetSetup_sim(n,nLevel)

    % -------- Get the true signal --------
    u=zeros(1400,1);
    u(91:210)= sin((0:119)/120*6*pi);
    u(251:550)= sin((0:299)/300*6*pi);
    u(651:1250)= sin((0:599)/600*6*pi);
    u(u<0) = 0;
    u = [zeros(100,1); u];
    x_exact = imresize(u,[n 1]);
    x_exact(x_exact<0.5)=0;
    x_exact(x_exact>=0.5)=1;
    output.x_exact = x_exact;
    
    % -------- Determine true impulse response function --------
    % Load h from text file
    h = textread('ImpulseResponse.txt');   % Put the name of the impulse response file here. It should be in txt format
    h = h(:,2);
    ht = [h; zeros(length(u)-length(h),1)];
    ht = imresize(ht,[n 1]);
    ht(ht<0) = 0; % h should be nonnegative
    ht = ht/sum(ht); % sum of h to be 1
    output.ht = ht;
    
    % -------- Define H matrix --------
    % Full matrix H
    Hfull = matrixH(ht,'zero');
    output.Hfull = Hfull;
    % Define matrix H by class
    H = matrixHfun(ht,'zero',1);
    output.H = H;
    
    % -------- Get delay and support of true h
    spt = ht>0.0001; 
    fndnz = find(spt==1); 
    delay = fndnz(1)-1; % delay
    output.delay = delay;
    
    supp = sum(spt); % number of entries in support
    output.supp = supp;
    
    % -------- true observation b
    b = H*x_exact;
    % adding noise to the output:
    rng(0)
    output.nLevel = nLevel;
    nsd = nLevel * norm(b) / sqrt(n);
    b = b + nsd*randn(size(b));
    output.b = b;
    
    % -------- Regularization matrix
    q11 = zeros(n,1); q11(1:2) = [1 -1]';
    q12 = zeros(n,1); q12(1) = 1;
    Q1 = toeplitz(q11,q12);
    
    q21 = zeros(n,1); q21(1:3)=[1 -2 1]';
    q22 = zeros(n,1); q22(1)=1;
    Q2 = toeplitz(q21,q22);      
    
    I = eye(n);
    
    % Choose regularization matrix Q
    output.Q = Q1;
    
end