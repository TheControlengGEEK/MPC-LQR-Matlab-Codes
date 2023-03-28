% All credits go to The Control eng GEEK.
% All rights reserved. Please mention The Control eng GEEK when 
% utilizing this code. Thanks.


clear all
close all
clc

% Process Model:
A = [ 0.4279  -0.2368   0.238    -0.07112  -0.04209
     -0.3167   0.4445   0.08599   0.0732    0.2731
      0.1261  -0.07832  0.2045    0.09331   0.2217
     -0.0244   0.1118   0.05136   0.7449   -0.1913
      0.05267  0.336    0.1149   -0.1822    0.5821];

B = [  0      -0.8637
      -0.1649  0.07736
       0      -1.214
       1.093  -1.114
       0      -0.006849];  
 
C = [ 1.533   0.3714    0       0    1.101
      0      -0.2256   -1.089   0    1.544]; 

% Cost function matrices:
q = eye(size(C,1))*1000; % weight on tracking error
s = eye(size(B,2))*0.01; % weight on control chnage
w = eye(size(C,1))*10;   % weight on slack variables
l = [0.1;0.2];

% For simulation: 
x0   = zeros(size(A,1),1);
r    = [1;3];               % tracking reference 
sh   = 20;                  % simiulation duration
umo  = zeros(size(B,2),1);  % control action at k = -1
x1 = []; x2 = []; x3 = []; x4 = []; x5 = [];
u1 = []; u2 = []; y1 = []; y2 = [];  

% For the MPC controller:
Np          = 5;  % Prediction Horizon
Nc          = Np; 
[PHI,Gamma] = MPC_Matrices(A,B,C,Np,Nc);

% Constraints:
su   = 0.5;   % on control move
du   = 0.05;  % on control change
yl   = 2;     % on controled output
posf = 10^5;  % infinity!
minf = -posf; % minus infinity!

Q = []; R = []; W = []; L = [];
  for i = 1:Np
      Q = blkdiag(Q,q);
      W = blkdiag(W,w);
  end
  for i = 1:Np
      R = [R;r];
      L = [L;l];
  end

Constraints_and_weights
  
% Setuping the MPC matrices:
  H  = blkdiag((Gamma'*Q*Gamma + PIE),W);
  Acon   = [gamma; -gamma]; 
  bcon   = [cmax ; -cmin]; 
  
  
  % running the MPC algorithm:
  for i = 1:sh
       f  = [(Gamma'*Q*(PHI*x0-R)-PSI*umo);L];
% MPC controller:
       U  = quadprog(H,f,Acon,bcon);                  
% The process under the action of MPC:
       x  = A*x0 + B*U(1:nu);
       y  = C*x0;
       
       x0 = x; umo = U(1:nu);
      dU_MIN = []; dU_MAX = [];
  for i = 1:Np
      dU_MIN = [dU_MIN;u_min]; dU_MAX = [dU_MAX;u_max];
  end
       dU_MIN(1:nu) = u_min + U(1:nu);
       dU_MAX(1:nu) = u_max + U(1:nu);
       cmin  = [[U_MIN;zeros((Np*ny),1)];     dU_MIN;  minf*ones(Np*ny,1);  (Y_MIN-PHI*x0);     zeros((Np*ny),1)];
       cmax  = [[U_MAX;ones((Np*ny),1)*posf]; dU_MAX; (Y_MAX-PHI*x0);       posf*ones(Np*ny,1);  posf*ones(Np*ny,1)];
       bcon  = [cmax ; -cmin];
       
% For ploting the outcomes:
       x1 = [x1;x(1)]; x2 = [x2;x(2)]; x3 = [x3;x(3)]; x4 = [x4;x(4)]; x5 = [x5;x(5)];
       u1 = [u1;U(1)]; u2 = [u2;U(2)];
       y1 = [y1;y(1)]; y2 = [y2;y(2)];
  end    

% Ploting the results: 
k = 0:sh;
figure(1)
subplot(3,1,1)
plot(k,[x0(1);x1],k,[x0(2);x2],k,[x0(3);x3],k,[x0(4);x4],k,[x0(5);x5])
legend('x1','x2','x3','x4','x5')

subplot(3,1,2)
stairs(k,[0;u1])
hold on
stairs(k,[0;u2])
hold off
legend('u1','u2')

subplot(3,1,3)
plot(k,[0;y1],k,[0;y2])
legend('y1','y2')

