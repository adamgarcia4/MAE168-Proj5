function [w, r, k] = timoshenkoElement(EA,EI, kGA, X, d, q)

%% Preprocessor Commands
clear all;
clc;

%% Find Transformation Matrix

X = [0 1 2 3 4 5];
d = [0 1 2 3 4 5]';
q = 1;
dx_ = X(4) - X(1);
dy_ = X(5) - X(2);
l_ = sqrt(dx_^2 + dy_^2);

c = dx_/l_;
s = dy_/l_;

T1 = [c -s 0; s c 0; 0 0 1];
T = [T1 zeros(3); zeros(3) T1];

%% Find K for element
syms x L;
N1 = 1- x/L;
N2 = x/L;

N1prime = diff(N1,x);
N2prime = diff(N2,x);

syms EI_ kGA_;

B = [0 N1prime 0 N2prime; -N1prime N1 -N2prime N2];
D = [EI_ 0;0 kGA_];

k_bef = transpose(B)*D*B;

k = int(k_bef,x,0,L);
%n=[N1 N2]
%int(n,x,0,l)

% syms x y;
% x+y
% f = x+y
% subs(f,x,3)
syms EA_;

k1 = EA_/L*[1 0 0 -1 0 0];
k2 = k(1,:);
k2 = [0 k2(1:2) 0 k2(3:4)];
k3 = k(2,:);
k3 = [0 k3(1:2) 0 k3(3:4)];
k4 = -k1;
k5 = k(3,:);
k5 = [0 k5(1:2) 0 k5(3:4)];
k6 = k(4,:);
k6 = [0 k6(1:2) 0 k6(3:4)];

% Complete Matrix (Symbolic)
k = [k1;k2;k3;k4;k5;k6];

%Numerical Substitutions
k = subs(k,EA_,5);
k = subs(k,EI_,3);
k = subs(k,L,1);
k = subs(k,kGA_,7)
%pretty(k)
%% Residual Force

fqprime = transpose([0 N1 0 0 N2 0]);
fqprime1 = int(q*fqprime,x,0,L);
fqprime1 = subs(fqprime1,L,l_);

r = k*(T'*d) - T * fqprime1


% % Distributed load vector
% fqPrime = [0 q*l_/2 q*l_^2/12 0 q*l_/2 -q*l_^2/12]';
% fq = T*fqPrime;
% 
% % Force in each element
% % rPrime = kPrime*(T'*d)-fqPrime;
% 
% r = k*d - fq;

%% Strain Energy
syms x L EI_ kGA_;
th1 = d(3);
th2 = d(6);
w1 = d(2);
w2 = d(5);
theta = th1*(1-x/L) + th2*(x/L);
thetaprime = diff(theta,x);
w_ = w1*(1-x/L)+w2*(x/L);
wprime = diff(w_,x);

w =  EI_*(thetaprime)^2 + kGA_*(wprime - theta)^2;
w1 = 1/2*int(w,x,0,L);
w2 = int(q*w_,x,0,L);

w1 = subs(w1,EI_,3);
w1 = subs(w1,L,1);
w1 = subs(w1,kGA_,7);

w2 = subs(w2,L,1);

w = w1-w2

end
