%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beamDiscussion.m
%
% Solves an EB frame using finite elements
%
% (c) 2015 MAE M168
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

% Example from section 4-7
L = 144;
X = [0   0  0;
     L   0  0;
     L+L/sqrt(2) -L/sqrt(2) 0]';
 
EQN = [0 1 0;
       0 2 0;
       0 3 0];
   
CNX = [1 2;
       2 3]';
   
E = [10*10^6; 10*10^6];
% I = [67.78; 67.78]; % Course reader states these values but rounds in solution
I = [68;68];
% A = [5.972; 5.972];
A = [6;6];
q = [-100;0];

kGA = [1; 1];

EI = E.*I;
EA = E.*A;

displacements = zeros(size(EQN));

[W, R, K] = beamAssembly(EA, EI, kGA, CNX, EQN, X, displacements,q);
W
R
K
 
Q = [0; 0; 0];

D = K\(Q-R);
D













