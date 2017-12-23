function [Oalp,Obeta]= Oracle_est(alpfun,betafun,u,sampleX,y,kC,kA,m,delta )
% construct oracle of component functions
% alpfun,betafun: true function for component functions
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline 
% sampleX,u,y are observations
% delta: ridge regression parameter

%oracle estimation for varying-coefficient functions
 T = length(y);
 Bt= rspline(u, u ,m, kC); 
 Oy1 = y - betafun(:,3);
 O1 = (kron(betafun(:,1), ones(1,size(Bt,2)))).*Bt;
 O2 = (kron(betafun(:,2), ones(1,size(Bt,2)))).*Bt;
 O3 = (kron(sampleX(:,4),ones(1,size(Bt,2)))).*Bt;
 O=[Bt O1 O2 O3];
 coeff1 = pinv(O'*O+delta*eye(size(O,2)))*O'*Oy1;
 B= blkdiag(Bt,Bt,Bt,Bt);
 alp = B *coeff1;
 alp = reshape(alp,T,4);
 Oalp = alp./repmat([1 mean(alp(:,2:end))],T,1);
 % oracle estimation for additive function 
 Oy2 = y - alpfun(:,1) - alpfun(:,5).*sampleX(:,4);
 Bx1 = rspline(sampleX(:,1),sampleX(:,1), m, kA);   
 Bx2 = rspline(sampleX(:,2),sampleX(:,2), m, kA);
 Bx3 = rspline(sampleX(:,3),sampleX(:,3), m, kA); 
 Bx1 = Bx1(:,2:end);
 Bx2 = Bx2(:,2:end);
 Bx3 = Bx3(:,2:end);
 Bx1 = Bx1 -repmat(mean(Bx1),T,1);
 Bx2 = Bx2 -repmat(mean(Bx2),T,1);
 Bx3 = Bx3 -repmat(mean(Bx3),T,1);
 OX1 = (kron(alpfun(:,2), ones(1,size(Bx1,2)))).*Bx1;
 OX2 = (kron(alpfun(:,3), ones(1,size(Bx1,2)))).*Bx2;
 OX=[OX1 OX2 Bx3];
 coeff2 =  pinv(OX'*OX+delta*eye(size(OX,2)))*OX'*Oy2;
 BX= blkdiag(Bx1,Bx2,Bx3);
 beta = BX*coeff2;
 Obeta = reshape(beta,T,3);
end

