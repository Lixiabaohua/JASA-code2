function [Oalp,Obeta]= Oracle_est(alpfun,betafun,u,sampleX,y,kC,kA,m,delta )
%construct oracle of component functions
% alpfun,betafun: true function for component functions
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline 
% sampleX,u,y are observations
% delta: ridge regression parameter

 T = length(y);
 Bt= rspline(u, u ,m, kC); 
 O1=(kron(betafun(:,1), ones(1,size(Bt,2)))).*Bt;
 O2=(kron(betafun(:,2), ones(1,size(Bt,2)))).*Bt;
 O=[Bt O1 O2];
 coeff1 = pinv(O'*O+delta*eye(size(O,2)))*O'*y;
 B= blkdiag(Bt,Bt,Bt);
 alp = B *coeff1;
 alp = reshape(alp,T,3);
 Oalp = alp./repmat([1 mean(alp(:,2:3))],T,1);
 % oracle estimation for additive function 
 Oy = y - alpfun(:,1);
 Bx1 = rspline(sampleX(:,1),sampleX(:,1), m, kA);   
 Bx2 = rspline(sampleX(:,2),sampleX(:,2), m, kA); 
 %Bx1 = Bx1(:,2:end);
 %Bx2 = Bx2(:,2:end);
 Bx1 = Bx1 - repmat(mean(Bx1),T,1);
 Bx2 = Bx2 - repmat(mean(Bx2),T,1);
 OX1 = (kron(alpfun(:,2), ones(1,size(Bx1,2)))).*Bx1;
 OX2 = (kron(alpfun(:,3), ones(1,size(Bx1,2)))).*Bx2;
 OX=[OX1 OX2];
 coeff2 =  pinv(OX'*OX+delta*eye(size(OX,2)))*OX'*Oy;
 BX= blkdiag(Bx1,Bx2);
 beta = BX*coeff2;
 Obeta = reshape(beta,T,2);
end

