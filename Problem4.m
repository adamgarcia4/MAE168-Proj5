%% Preprocessor Commands

clear all;
clc;

%% Parameters
theta_total=90;
n_el = 10;
n_nodes = n_el+1;
radius = 1;
EA_const = 1e9;
EI_const = 1e9;

%% Prepadding for Efficiency
% X contains all node element positions in (x,y,z) format
EA = ones(n_el,1);
EI = ones(n_el,1);
X = zeros(n_nodes,3);
CNX = zeros(2,n_el);

%% Generate EA/EI arrays
EA = EA_const*EA; %Make a vector here and below.
EI = EI_const*EI;

%% Generate CNX

for iter = 1:n_el
    CNX(1,iter) = iter;
    CNX(2,iter) = iter+1;
end
CNX;

%% Generate EQN

counter = 1;
EQN = zeros(3,n_nodes);
for i = 1:n_nodes
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


%% Generate X
ThetaArray = linspace(0,theta_total*pi/180,n_nodes);

for iter = 1:n_nodes
    X(iter, 1) = -radius*cos(ThetaArray(iter));
    X(iter, 2) = radius*sin(ThetaArray(iter));
    X(iter, 3) = pi/2-ThetaArray(iter);
end
X;

% displ = zeros(size(EQN);
% [W, R, K] = beamAssembly(EA, EI, CNX, EQN, X, displ, q);
% 
% Q = [0;0;0];
% D = K\(Q-R);
% D
% 
% 




