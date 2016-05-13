%% Preprocessor Commands

clear all;
close all;
clc;
fignum = 1;
exist([pwd,'\plots'],'dir');% check whether a folder named 'plots' exists
% in the working directory; returns the value 7
% if such a directory exists, 0 otherwise
if ans == 0; % if the folder 'plots' does not exist in the working directory
mkdir('plots') % create a folder in which to store images
end

%% Inputs

% (EA,EI, kGA, X, d, q)

EA = 60*10^6;
EI = 680*10^6;
kGA = 11603*6*1;
number = 10;
% EI/kGAL^2 vary this dimensionless number.
x = ones(number,2);
Q = [0;-100; 0];
q = 0;

EQN = [0 1;
       0 2;
       0 3];
CNX = [1 2]';

p1 = zeros(2,number);


for iter = 1:number
    if iter == 1
        L = sqrt(EI/(kGA*iter))
    elseif iter == 100
        L = sqrt(EI/(kGA*iter))
    else
        L = sqrt(EI/(kGA*iter));
    end    
    X = [0 0 0
         L 0 0]';
     dt = zeros(6,1);
     [Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,dt,q);
     dt = Kt\(Q);
     
     
     p1(1,iter) = iter;
     p1(2,iter) = dt(2);
end

p2 = zeros(2,number);
for iter = 1:number
    
    L = sqrt(EI/(kGA*iter));
    X = [0 0 0
         L 0 0]';
    d = zeros(6,1);
    [W, R, K]= beamAssembly( EA, EI, CNX, EQN, X, d, q);
    d = K\Q;
    
     p2(1,iter) = iter;
     p2(2,iter) = d(2);
    
    
end

figure(1)
semilogx(fliplr(p1(1,:)),p1(2,:)+.04,'LineWidth',2)
hold on
semilogx(fliplr(p2(1,:)),p2(2,:),'r','LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection on Cantilevered Beam 1-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_1 ','Tip Deflection on Cantilevered Beam 1-element','.png']);

figure(2)
diff = -fliplr(p1(2,:)) + fliplr(p2(2,:))
semilogx((p1(1,:)),diff,'LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection Difference','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection Difference on Cantilevered Beam 1-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_2 ','Tip Deflection Difference on Cantilevered Beam 1-element','.png']);
