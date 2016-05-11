%% Preprocessor Commands
clear all;
clc;

%% Find Transformation Matrix

x = [0 1 2 3 4 5];
dx_ = x(4) - x(1);
dy_ = x(5) - x(2);
l_ = sqrt(dx_^2 + dy_^2);

c = dx_/l_;
s = dy_/l_;

T1 = [c -s 0; s c 0; 0 0 1];
T = [T1 zeros(3); zeros(3) T1];

%% Find K for element
syms x l;
N1 = 1- x/l;
N2 = x/l;

N1prime = diff(N1,x);
N2prime = diff(N2,x);

syms EI kGA;

B = [0 N1prime 0 N2prime; -N1prime N1 -N2prime N2];
D = [EI 0;0 kGA];
%pretty(B)
%pretty(D)

k_bef = transpose(B)*D*B;

%pretty(k_bef);

k = int(k_bef,x,0,l);
%n=[N1 N2]
%int(n,x,0,l)

% syms x y;
% x+y
% f = x+y
% subs(f,x,3)
%% hi
syms EA;

k1 = EA/l*[1 0 0 -1 0 0];
k2 = k(1,:);
k2 = [0 k2(1:2) 0 k2(3:4)];
k3 = k(2,:);
k3 = [0 k3(1:2) 0 k3(3:4)];
k4 = -k1;
k5 = k(3,:);
k5 = [0 k5(1:2) 0 k5(3:4)];
k6 = k(4,:);
k6 = [0 k6(1:2) 0 k6(3:4)];
k = [k1;k2;k3;k4;k5;k6]
pretty(k)

%% Residual Force

r = k*d - q
