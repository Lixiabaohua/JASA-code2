function [mis_alp,mis_beta,Vcoeff,Acoeff, para ] = mis_Pest(sampleX,y,Bt,Bx1,Bx2,delta )
%misspecified additive model and varying-coefficient model
%estiamtion of misespecified varying-coefficient model
% Bx1,Bx2,Bt,b1,b2:B-spline basis matrix
% delta: ridge regression parameters
% Vcoeff: coefficients for varying-coefficient functions
% Acoeff: coefficients for additive functions
T=length(y);
V1 = (kron(sampleX(:,1), ones(1,size(Bt,2)))).*Bt;
V2=  (kron(sampleX(:,2), ones(1,size(Bt,2)))).*Bt;
V=[Bt V1 V2];
Vcoeff = pinv(V'*V+delta*eye(size(V,2)))*V'*y;
B = blkdiag(Bt,Bt, Bt);
alpvec = B*Vcoeff;
mis_alp = reshape(alpvec,T,3);
%estiamtion of misespecified additive model
%\alpha_{1}=\alpha_{2}=1,\alpha_{0} is unknown constant
Bx = [ones(T,1),Bx1 Bx2];
Acoeff =  pinv(Bx'*Bx+delta*eye(size(Bx,2)))*Bx'*y;
BX = blkdiag(Bx1,Bx2);
betavec = BX*Acoeff(2:end);
ebeta = reshape(betavec,T,2);
meanvec = mean(ebeta);
ebeta = ebeta - repmat(meanvec,T,1);
const =  Acoeff(1) + meanvec(1)  +meanvec(2) ;
mis_beta = [const*ones(T,1) ebeta];
para = [meanvec const];
end

