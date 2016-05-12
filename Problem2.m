%% Preprocessor Commands

clear all;
clc;

%% Inputs

% (EA,EI, kGA, X, d, q)

EA = 1;
EI = 1;
kGA = 1;
number = 100;
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
Q = [0;-100; 0];
q = 0;

EQN = [0 1;
       0 2;
       0 3];

CNX = [1 2]';
p = zeros(2,100);
for iter = 1:100
    L = sqrt(EI/(kGA*iter));
    X = [0 0 0
         L 0 0];
     dt = zeros(6,1);
     [Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,d,q);
     dt = Kt\(Q);
     p(1,iter) = iter;
     p(2,iter) = dt(2);
end

plot(p(1,:),p(2,:))


%% Computations