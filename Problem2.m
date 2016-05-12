%% Preprocessor Commands

clear all;
clc;

%% Inputs

% (EA,EI, kGA, X, d, q)

EA = 1;
EI = 1;
kGA = 1;
number = 100
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
for iter = 1:100
    L = sqrt(EI/(kGA*iter));
    X = [0 0 0
         L 0 0];
     x(iter,1) = iter*x(iter,1);
     x(iter,2) = iter*x(iter,2);
     
     
end

plot(x(:,1),x(:,2))


%% Computations