%% Preprocessor Commands

clear all;
close all;
clc;

%% Inputs

% (EA,EI, kGA, X, d, q)

EA = 60*10^6;
EI = 680*10^6;
kGA = 11603*6*1;
number = 10;
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
Q = [0;-100; 0];
q = 0;

EQN = [0 1;
       0 2;
       0 3];
CNX = [1 2]';

p1 = zeros(2,number);


for iter = 1:number
    if iter == 1
        L = sqrt(EI/(kGA*iter))
    elseif iter == 100
        L = sqrt(EI/(kGA*iter))
    else
        L = sqrt(EI/(kGA*iter));
    end
    X = [0 0 0
         L 0 0]';
     dt = zeros(6,1);
     [Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,dt,q);
     dt = Kt\(Q);
     
     
     p1(1,iter) = iter;
     p1(2,iter) = dt(2);
end

p2 = zeros(2,number);
for iter = 1:number
    
    L = sqrt(EI/(kGA*iter));
    X = [0 0 0
         L 0 0]';
    d = zeros(6,1);
    [W, R, K]= beamAssembly( EA, EI, CNX, EQN, X, d, q);
    d = K\Q;
    
     p2(1,iter) = iter;
     p2(2,iter) = d(2);
    
    
end


figure(1)
semilogx(p1(1,:),p1(2,:))
hold on
semilogx(p2(1,:),p2(2,:))




%% Computations