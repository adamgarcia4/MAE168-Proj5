%% Preprocessor Commands

clear all;
clc;
exist([pwd,'\plots'],'dir');% check whether a folder named 'plots' exists
% in the working directory; returns the value 7
% if such a directory exists, 0 otherwise
if ans == 0; % if the folder 'plots' does not exist in the working directory
mkdir('plots') % create a folder in which to store images
end

%% Parameters
theta_total=90;
n_el = 10;
number = n_el;
n_nodes = n_el+1;
radius = 1;
EA_const = 1e9;
EI_const = 1e9;

%% Prepadding for Efficiency
% X contains all node element positions in (x,y,z) format
EA = 60*10^6*ones(1,10)';
EI = 680*10^6*ones(1,10)';
kGA = 11603*6*1*ones(1,10)';
q = zeros(10,1);

Q = zeros(30,1);
Q(29,1) = -100;

X = zeros(n_nodes,3);
CNX = zeros(2,n_el);



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
X = transpose(X);

%% Running loops
p1 = zeros(2,number);
for iter = 1:number
    if iter == 1 || iter==number
        l = sqrt(EI(1)/(kGA(1)*iter))
    else
        l = sqrt(EI(1)/(kGA(1)*iter));
    end
    
     dt = zeros(size(EQN));
     [Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,dt,q);
     dt = Kt\(Q);
     
     
     p1(1,iter) = iter;
     p1(2,iter) = dt(29);
end

p2 = zeros(2,number);
for iter = 1:number
    
    l = sqrt(EI(1)/(kGA(1)*iter));
    
    d = zeros(size(EQN));
    [W, R, K]= beamAssembly( EA, EI, CNX, EQN, X, d, q);
    d = K\Q;
    
     p2(1,iter) = iter;
     p2(2,iter) = d(29);
    
    
end

%% Plotting
figure(5)
semilogx(fliplr(p1(1,:)),p1(2,:)+.04,'LineWidth',2)
hold on
semilogx(fliplr(p2(1,:)),p2(2,:),'r','LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection on Curved Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_5 ','Tip Deflection on Curved Cantilevered Beam 10-element','.png']);

figure(6)
diff = -fliplr(p1(2,:)) + fliplr(p2(2,:))
semilogx((p1(1,:)),diff,'LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection Difference','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection Difference on Curved Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_6 ','Tip Deflection Difference on Curved Cantilevered Beam 10-element','.png']);


