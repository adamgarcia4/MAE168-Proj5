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
radius = 1;

%% Prepadding for Efficiency
% X contains all node element positions in (x,y,z) format

n_el_arr = 2:10:50;


%% Start For Loop
for iter = 1:size(n_el_arr,2)
n_el = n_el_arr(iter); 
number = n_el;
n_nodes = n_el +1;
EA = 60*10^6*ones(1,n_el)';
EI = 680*10^6*ones(1,n_el)';
kGA = 11603*6*1*ones(1,n_el)';


q = zeros(n_el,1);
Q = zeros(3*n_el,1);
Q(3*n_el-1,1) = -100;
X = zeros(n_nodes,3);
CNX = zeros(2,n_el);

%% Generate CNX
for i = 1:n_el
    CNX(1,i) = i;
    CNX(2,i) = i+1;
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

for ind = 1:n_nodes
    X(ind, 1) = -radius*cos(ThetaArray(ind));
    X(ind, 2) = radius*sin(ThetaArray(ind));
    X(ind, 3) = pi/2-ThetaArray(ind);
end
X = transpose(X);

%% Run Loops

l = sqrt((X(1,1)-X(2,1))^2+(X(1,2)-X(2,2))^2);
d = zeros(size(EQN));
[Wt, Rt, Kt] = beamTimoshenkoAssembly(EA,EI,kGA,CNX, EQN,X,d,q);
d = Kt\(Q);
p1(1,iter) = l; %plotting against element length
p1(2,iter) = d(3*n_el-1); %Tip Transverse Deflection
p1(3,iter) = d(3*n_el); %Tip Rotation
end

%% Running loops


% p2 = zeros(2,number);
% for iter = 1:number
%     
%     l = sqrt(EI(1)/(kGA(1)*iter));
%     
%     d = zeros(size(EQN));
%     [W, R, K]= beamAssembly( EA, EI, CNX, EQN, X, d, q);
%     d = K\Q;
%     
%      p2(1,iter) = iter;
%      p2(2,iter) = d(29);
%     
%     
% end

%% Plotting
figure(5)
plot(fliplr(p1(1,:)),p1(2,:),'LineWidth',2)
xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
ylabel('Deflection','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
title('Tip Deflection on Curved Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
%grid on;
set(gca, 'FontSize', 16);
print('-dpng',[pwd,'\plots\Figure_5 ','Tip Deflection on Curved Cantilevered Beam 10-element','.png']);

% figure(6)
% diff = -fliplr(p1(2,:)) + fliplr(p2(2,:))
% semilogx((p1(1,:)),diff,'LineWidth',2)
% xlabel('log(L/h)','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
% ylabel('Deflection Difference','FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Times New Roman')
% title('Tip Deflection Difference on Curved Cantilevered Beam 10-element', 'FontSize', 20, 'FontName', 'Times New Roman');
% %grid on;
% set(gca, 'FontSize', 16);
% print('-dpng',[pwd,'\plots\Figure_6 ','Tip Deflection Difference on Curved Cantilevered Beam 10-element','.png']);
% 
% 
