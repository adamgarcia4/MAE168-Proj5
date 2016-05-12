%% Preprocessor Commands
clear; close all; clc;

%% Inputs
%(EA, EI, kGA, CNX, EQN, X, displacements,q);
EA = 60*10^6;
EI = 680*10^6;
kGA = 11603*6*1;
number = 1000;
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
Q = [0;-100; 0];
q = 0;
numElements = 10;
numNodes = 11;

%% Generate EQN
EQN = zeros(3,numNodes);
counter = 1;
for i = 1:numNodes
   if i ==1
       EQN(:,1)=0;
   else
       EQN(1,i)=counter;
       counter = counter+1;
       EQN(2,i)=counter;
       counter = counter+1;
       EQN(3,i)=counter;
       counter = counter+1;
   end
end

%% Generate CNX

CNX = zeros(numElements,2);
counter = 1;
for i = 1:numElements
   CNX(i,1) = counter;
   CNX(i,2) = counter + 1;
   counter = counter + 1;    
end

%% Generate Q

Q = zeros(1,numNodes);
Q(1,11) = -100;

%% Generate X

X = zeros(numNodes,3);
l = 
for i = 1:numNodes
   X(i,1) =  
    
    
end

%% Generate EQN and CNX

%(X, displacements);
%(EA, EI, kGA, CNX, EQN, X, displacements,q);














% 
% % Example from section 4-7
% L = 144;
% X = [0   0  0;
%      L   0  0;
%      L+L/sqrt(2) -L/sqrt(2) 0]';
%  
% EQN = [0 1 0;
%        0 2 0;
%        0 3 0];
%    
% CNX = [1 2;
%        2 3]';
%    
% E = [10*10^6; 10*10^6];
% % I = [67.78; 67.78]; % Course reader states these values but rounds in solution
% I = [68;68];
% % A = [5.972; 5.972];
% A = [6;6];
% q = [-100;0];
% 
% kGA = [1; 1];
% 
% EI = E.*I;
% EA = E.*A;
% 
% displacements = zeros(size(EQN));
% 
% [W, R, K] = beamTimoshenkoAssembly(EA, EI, kGA, CNX, EQN, X, displacements,q);
% % W
% R
% K
%  
% Q = [0; 0; 0];
% 
% D = K\(Q-R);
% % D
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
