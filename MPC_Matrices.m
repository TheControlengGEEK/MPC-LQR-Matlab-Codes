function [PHI,Gamma] = MPC_Matrices(A,B,C,Nc,Np)

[~,ru] = size(B);
[py,~] = size(C);

% Creating capital PHI:
hy  = C;
PHI = C*A;

for k = 2:Np
hy(py*(k-1)+1:py*k,:) = hy(py*(k-1)+1-py:py*(k-1),:)*A;
PHI(py*(k-1)+1:py*k,:) = PHI(py*(k-1)+1-py:py*(k-1),:)*A;
end

% Creating capital Gamma:
v1                   = hy*B; 
[~,mv1]              = size(v1);
Gamma                = zeros(py*Np,ru*Nc);
Gamma(1:py*Np,1:mv1) = v1;

for i = 2:Nc              
Gamma(:,mv1*(i-1)+1:mv1*i) = [zeros((i-1)*py,mv1);v1(1:(Np-i+1)*py,1:mv1)];   
end

