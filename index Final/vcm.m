function [ halp,res,sig ] =vcm(Sx,u,Sy,m,k,delta)
% varying-coefficient model estimation
% Sx and Sy: for sample of covariates and response
% m: order of B-spline
% k: knot number

T = length(Sy);
Bt= rspline(u, u ,m, k); 
D1=(kron(Sx(:,1), ones(1,size(Bt,2)))).*Bt;
D2=(kron(Sx(:,2), ones(1,size(Bt,2)))).*Bt;
D=[Bt D1 D2 ];
coeff1 = pinv(D'*D+delta*eye(size(D,2)))*D'*Sy;
B= blkdiag(Bt,Bt,Bt);
halp = B*  coeff1;
halp = reshape(halp,T,3);
fit = halp(:,1) + halp(:,2) .* Sx(:,1) +  halp(:,3) .* Sx(:,2);
res = Sy - fit;
p = 3 * (m + k);
sig = sum(res.^2)/(T-p);
end

