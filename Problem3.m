%% Preprocessor Commands
clear; close all; clc;

exist([pwd,'\plots'],'dir');% check whether a folder named 'plots' exists
% in the working directory; returns the value 7
% if such a directory exists, 0 otherwise
if ans == 0; % if the folder 'plots' does not exist in the working directory
mkdir('plots') % create a folder in which to store images
end

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

figure(3)
semilogx(fliplr(p1(1,:)),p1(2,:)+.04,'LineWidth',2)
hold on
semilogx(fliplr(p2(1,:)),p2(2,:),'r','LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection on Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_3 ','Tip Deflection on Cantilevered Beam 10-element','.png']);

figure(4)
diff = -fliplr(p1(2,:)) + fliplr(p2(2,:))
semilogx((p1(1,:)),diff,'LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection Difference','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection Difference on Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_4 ','Tip Deflection Difference on Cantilevered Beam 10-element','.png']);
