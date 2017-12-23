function [Valp,Abeta]= Mis_est( u,sampleX,y,kC,kA,m,delta )
%estiamtion of misespecified varying-coefficient model
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline in Step II and III estimation
% sampleX, u, y:sample
% delta: redge regression parameter
T=length(y);
Bt= rspline(u, u ,m, kC); 
V1 = (kron(sampleX(:,1), ones(1,size(Bt,2)))).*Bt;
V2=  (kron(sampleX(:,2), ones(1,size(Bt,2)))).*Bt;
V=[Bt V1 V2];
Vcoeff = pinv(V'*V+delta*eye(size(V,2)))*V'*y;
B= blkdiag(Bt,Bt,Bt);
alp = B *Vcoeff;
Valp = reshape(alp,T,3);
%estiamtion of misespecified additive model
%\alpha_{1}=\alpha_{2}=1,\alpha_{0} is unknown constant
Bx1 = rspline(sampleX(:,1),sampleX(:,1), m, kA);   
Bx2 = rspline(sampleX(:,2),sampleX(:,2), m, kA);
Bx = [ones(T,1),Bx1 Bx2];
coeff =  pinv(Bx'*Bx+delta*eye(size(Bx,2)))*Bx'*y;
BX= blkdiag(Bx1,Bx2);
beta = BX*coeff(2:end);
beta = reshape(beta,T,2);
center = mean(beta);
const =  coeff(1) +center(1)+center(2);
Abeta =[ repmat(const,T,1) beta];
end

