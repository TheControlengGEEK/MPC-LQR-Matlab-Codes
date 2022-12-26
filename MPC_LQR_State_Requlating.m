% All credits go to The Control eng GEEK.
% All rights reserved. Please mention The Control eng GEEK when 
% utilizing this code. Thanks.


clear all
close all
clc

% Process Model:
A = [-2.3122   -0.0093   -2.1478    0.1603    1.3486
     -0.2226   -0.2953    2.1721    2.1696   -3.3464
      0.2684   -1.0571   -1.1903    0.0804    2.0020
      2.0505    2.0088   -0.1674   -1.2369   -1.1968
     -0.2121    1.4738   -0.0478    0.0920   -1.2225];

B = [ 0.0000   -0.5900   
      2.2294   -0.2781  
      0.3376    0.4227  
      1.0001   -1.6702  
      0.0000    0.4716];  
 
C = [ 0.0000    0.2571   -1.3218    0.0000    0.9111
     -0.6509   -0.9444    0.0000    0.0000    0.5946]; 

D = zeros(size(C,1),size(B,2)); 

% Cost function matrices:
q = eye(size(A))*100;
s = eye(size(B,2))*0.01;

% Checking controllability and observability:
unco = length(A) - rank(ctrb(A,B));  
unob = length(A) - rank(obsv(A,C));

% Computing the terminal stablizing weight P:
nc  = size(A,1);   nu  = size(B,2);
qsc = sqrtm(q);    rsc = sqrtm(s);   
Z1  = zeros(nc,nc); Z2 = zeros(nc,nu); I1  = eye(nc); I2 = eye(nu);

cvx_begin sdp
variable E(nc,nc) symmetric; 
variable L(nu,nc);
E >= 0; 

[       E    (A*E+B*L)'   E*qsc    L'*rsc;            
    A*E+B*L      E         Z1        Z2;
     qsc*E       Z1'       I1        Z2;
     rsc*L      Z2'        Z2'       I2] >= 0;
cvx_end

P = inv(E);
K = L*P;
Aclosed = (A + B*K); % because we assume u = Kx in the LQR contoller.

% For simulation: 
x0 = [0.2925; 0.2177;0.3019;0.5954; 0.5637];
xp = x0;
u1 = []; u2 = []; 
x1 = []; x2 = []; x3 = []; x4 = []; x5 = []; 
sh = 20; % simiulation duration

% Desiging the Kalman filter to estimate the states from the inputs and outputs:
Gp            = eye(size(A)); % assume all of the states subjected to noise 
sys           = ss(A,[B Gp],C,[D zeros(2,5)],-1);
p_noise       = normrnd(0,0.00008,sh,size(A,1)); 
m_noise       = normrnd(0,0.00005,sh,2); 
Qw            = cov(p_noise); 
Rv            = cov(m_noise);
[kalmf,Lc,Pc] = kalman(sys,Qw,Rv); 
eigKMLc       = eig(A - Lc*C);


% For the MPC controller:
Np      = 3;  
[Fx,Mx] = MPC_Matrices(A,B,eye(size(A)),Np,Np);
su      = 0.5;
uv_min  = -su*ones(size(B,2),1); uv_max = su*ones(size(B,2),1); 
UV_MIN  = []; UV_MAX = [];  
  for i = 1:Np
      UV_MIN  = [UV_MIN ; uv_min];
      UV_MAX  = [UV_MAX ; uv_max];
  end
Iin    = eye(Np*size(B,2)); 
Acon   = [Iin    ; -Iin]; 
bcon   = [UV_MAX ; -UV_MIN];
Q = []; S = [];
  for i = 1:Np-1
      Q = blkdiag(Q,q);
  end
  for i = 1:Np
      S = blkdiag(S,s);
  end
  Qp   = blkdiag(Q,P); % Combine the LQR with MPC (terminal weight)
  
% Running the simulation along with the MPC-LQR conrtroller & Kalman filter:
       H  = (Mx'*Qp*Mx + S);    
  for i = 1:sh
       f  = (Mx'*Qp*Fx*xp);
% MPC controller:
       U  = quadprog(H,f,Acon,bcon);                  
% The process under the action of MPC and noises:
       x  = A*xp + B*U(1:size(B,2)) + Gp*p_noise(i,:)';
       y  = C*xp + m_noise(i,:)';
% Estimating the states from the inputs and outputs:
       xe = A*xp + B*U(1:size(B,2)) + Lc*(y - C*xp);   % Kalman filter 
       xp = xe;    % feedback the setimated states
% For ploting the outcomes:
       x1 = [x1;x(1)]; x2 = [x2;x(2)]; x3 = [x3;x(3)];
       x4 = [x4;x(4)]; x5 = [x5;x(5)];
       u1 = [u1;U(1)]; u2 = [u2;U(2)];
  end    

% Ploting the results: 
k = 0:sh;
figure(2)
subplot(2,1,1)
plot(k,[x0(1);x1],k,[x0(2);x2],k,[x0(3);x3],k,[x0(4);x4],k,[x0(5);x5])
legend('x1','x2','x3','x4','x5')

subplot(2,1,2)
stairs(k,[0;u1])
hold on
stairs(k,[0;u2])
hold off
legend('u1','u2')

