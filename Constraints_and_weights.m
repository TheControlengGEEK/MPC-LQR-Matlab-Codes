nu = size(B,2);
ny = size(C,1);

u_min  = -su*ones(nu,1); 
u_max  =  su*ones(nu,1); 
du_min = -du*ones(nu,1); 
du_max =  du*ones(nu,1); 
y_min =   yl*ones(ny,1); 
y_max =  -yl*ones(ny,1);

U_MIN  = []; U_MAX  = []; dU_MIN = []; dU_MAX = []; Y_MIN = []; Y_MAX = [];
  for i = 1:Np
      U_MIN  = [U_MIN;u_min];  U_MAX  = [U_MAX;u_max];
      dU_MIN = [dU_MIN;u_min]; dU_MAX = [dU_MAX;u_max];
      Y_MIN  = [Y_MIN;y_min];  Y_MAX  = [Y_MAX;y_max];
  end
   
% for constraints on U:
Ibar = blkdiag(eye(Np*nu),eye(Np*nu));

% for constraints on delta U:
Omega = zeros((Np-1)*nu,Np*nu); % forming capital omega:
O = [-eye(nu,nu) eye(nu,nu)];
for k = 1:Np-1
 Omega((k-1)*nu+1:k*nu,(k-1)*nu+1:(k+1)*nu) = O;
end
Omega    = [eye(nu,nu) zeros(nu,size(Omega,2)-nu);Omega]; 
Omegabar = [Omega zeros(Np*ny,Np*ny)];

% for constriants on the controlled output z:
Gammabar1 = [Gamma -eye(Np*ny)];
Gammabar2 = [Gamma  eye(Np*ny)];

% for constriants on the slack variables:
Gammabar3 = [zeros(Np*nu,Np*ny)  eye(Np*ny)];

% combinig all constriants:
gamma = [Ibar;Omegabar;Gammabar1;Gammabar2;Gammabar3];
cmin  = [[U_MIN;zeros((Np*ny),1)];     dU_MIN;  minf*ones(Np*ny,1);  (Y_MIN-PHI*x0);     zeros((Np*ny),1)];
cmax  = [[U_MAX;ones((Np*ny),1)*posf]; dU_MAX; (Y_MAX-PHI*x0);       posf*ones(Np*ny,1);  posf*ones(Np*ny,1)];

% for Pi:
PIE = zeros(Np*nu,Np*nu);
if Np == 1
   PIE = s ;
else
   k = 0;
   PIE(1:nu,1:nu) = 2*s;
   PIE(1+nu:nu+nu,1:nu) = -s;

   for k = 1:Np-2
       ku = k*nu ;
       PIE(ku-nu+1:ku,ku+1:ku+nu) = -s;
       PIE(ku+1:ku+nu,ku+1:ku+nu) = 2*s;
       PIE(ku+nu+1:ku+2*nu,ku+1:ku+nu) = -s;
   end

   k = Np-1;
   ku = k*nu;
   PIE(ku-nu+1:ku,ku+1:ku+nu) = -s;
   PIE(ku+1:ku+nu,ku+1:ku+nu) = s;
end

% for Psi:
PSI = [-s;zeros((Np-1)*nu,nu)];



