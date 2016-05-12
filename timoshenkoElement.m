function [w, r, k] = timoshenkoElement(EA,EI, kGA, X, d, q)

%% Preprocessor Commands
clear all;
clc;

%% Find Transformation Matrix

X = [0 1 2 3 4 5];
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

syms ei kGA;

B = [0 N1prime 0 N2prime; -N1prime N1 -N2prime N2];
D = [ei 0;0 kGA];

k_bef = transpose(B)*D*B;

k = int(k_bef,x,0,L);
%n=[N1 N2]
%int(n,x,0,l)

% syms x y;
% x+y
% f = x+y
% subs(f,x,3)
syms ea;

k1 = ea/L*[1 0 0 -1 0 0];
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
k = subs(k,ea,5);
k = subs(k,ei,3);
k = subs(k,L,1);
k = subs(k,kGA,7);
%pretty(k)
%% Residual Force

% Distributed load vector
fqPrime = [0 q*l/2 q*l^2/12 0 q*l/2 -q*l^2/12]';
fq = T*fqPrime;

% Force in each element
% rPrime = kPrime*(T'*d)-fqPrime;

r = k*d - fq;

%% Strain Energy
syms x L ei kga;
th1 = X(3);
th2 = X(6);
w1 = X(2);
w2 = X(5);
theta = th1*(1-x/L) + th2*(x/L);
thetaprime = diff(theta,x);
w_ = w1*(1-x/L)+w2*(x/L);
wprime = diff(w_,x)

w =  ei*(thetaprime)^2 + kga*(wprime - theta)^2;
w = 1/2*int(pi,x,0,L)
end
