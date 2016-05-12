%% Preprocessor Commands
clear; close all; clc;

%% Inputs
%(EA, EI, kGA, CNX, EQN, X, displacements,q);
EA = 60*10^6*ones(1,10)';
EI = 680*10^6*ones(1,10)';
kGA = 11603*6*1*ones(1,10)';
number = 10;
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
q = zeros(10,1);
numElements = 10;
numNodes = numElements+1;

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

CNX = zeros(2,numElements);
counter = 1;
for i = 1:numElements
   CNX(1,i) = counter;
   CNX(2,i) = counter + 1;
   counter = counter + 1;    
end
% CNX = transpose(CNX);

%% Generate Q

Q = zeros(30,1);
Q(29,1) = -100;

%% Generate X

syms L Lcount X
Lcount=0;
for i=1:numNodes 
        for j=1:3
            if j==2 || j ==3
                X(j,i)=0;
            else
                X(j,i)= Lcount;
                Lcount = Lcount +L/numElements;
            end
        end
end
Xunm = X;
clear X;

%% Run through dimensionless params

p1 = zeros(2,number);
for iter = 1:number
    if iter == 1 || iter==number
        l = sqrt(EI(1)/(kGA(1)*iter))
    else
        l = sqrt(EI(1)/(kGA(1)*iter));
    end
    
    X = subs(Xunm,L,l);
    X = double(X);
     dt = zeros(size(EQN));
     [Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,dt,q);
     dt = Kt\(Q);
     
     
     p1(1,iter) = iter;
     p1(2,iter) = dt(29);
end

p2 = zeros(2,number);
for iter = 1:number
    
    l = sqrt(EI(1)/(kGA(1)*iter));
    X = subs(Xunm,L,l);
    X = double(X);
    
    d = zeros(size(EQN));
    [W, R, K]= beamAssembly( EA, EI, CNX, EQN, X, d, q);
    d = K\Q;
    
     p2(1,iter) = iter;
     p2(2,iter) = d(29);
    
    
end



figure(1)
semilogx(p1(1,:),p1(2,:))
hold on
semilogx(p2(1,:),p2(2,:))











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
