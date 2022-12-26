function [Fy,My] = MPC_Matrices(A,Bu,Cy,Nc,Np)

% Create Fy & My
[n,n] = size(A);
[n,ru] = size(Bu);
[py,n] = size(Cy);

hy = Cy;
Fy = Cy*A;
for k = 2:Np
hy(py*(k-1)+1:py*k,:) = hy(py*(k-1)+1-py:py*(k-1),:)*A;
Fy(py*(k-1)+1:py*k,:) = Fy(py*(k-1)+1-py:py*(k-1),:)*A;
end

v1 = hy*Bu; 
[nv1,mv1] = size(v1);
My = zeros(py*Np,ru*Nc);
My(1:py*Np,1:mv1) = v1;

for i = 2:Nc              
My(:,mv1*(i-1)+1:mv1*i) = [zeros((i-1)*py,mv1);v1(1:(Np-i+1)*py,1:mv1)];   
end

