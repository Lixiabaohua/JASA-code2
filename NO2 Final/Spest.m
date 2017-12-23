function [halp,hbeta,fit,res,sig] =Spest( Inibeta,k1, k2,m, sampleX, u, y,delta)
%Step II and Step III estimation (using same order of B-spline
% k1: knots for varying-coefficient function
%k2:  knots for additive function
%m: order of B-spline
T = length(y);
% Step II B-spline to estimate varying-coefficient functions
 Bt= rspline(u, u ,m, k1); 
 D1=(kron(Inibeta(:,1), ones(1,size(Bt,2)))).*Bt;
 D2=(kron(Inibeta(:,2), ones(1,size(Bt,2)))).*Bt;
 D3=(kron(Inibeta(:,3), ones(1,size(Bt,2)))).*Bt;
 D=[Bt D1 D2 D3];
 coeff1    = pinv(D'*D+delta*eye(size(D,2)))*D'*y;
 %estimation fot varying-coefficient function
 B= blkdiag(Bt,Bt,Bt,Bt);
 halp = B*  coeff1;
 halp = reshape(halp,T,4);
 %standardized
 mynorm = mean(halp(:,2:4));
 halp = halp./repmat([1 mynorm],T,1);
% Step III B-spline  to estimate additive functions
 Bx1 = rspline(sampleX(:,1),sampleX(:,1), m,k2);   
 Bx2 = rspline(sampleX(:,2),sampleX(:,2), m,k2);   
 Bx3 = rspline(sampleX(:,3),sampleX(:,3), m,k2);     
 pseuy =y - halp(:,1) ;
 S1 =  (kron(halp(:,2), ones(1,size(Bx1,2)))).*Bx1;
 S2 =  (kron(halp(:,3), ones(1,size(Bx2,2)))).*Bx2;
 S3 =  (kron(halp(:,4), ones(1,size(Bx3,2)))).*Bx3;
 S   = [S1 S2 S3];
 coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S'*pseuy;
 %estimation for additive functions
 BX= blkdiag(Bx1,Bx2,Bx3);
 hbeta = BX *coeff2;
 hbeta = reshape(hbeta,T,3);
 %centering
 c =mean(hbeta);
 hbeta = hbeta-repmat(c,T,1);
 %updated estimation for alpha0
 up0 = halp(:,1)+  halp(:,2)*c(1) +halp(:,3) *c(2)+halp(:,4 ) * c(3);
 fit = up0 + halp(:,2).* hbeta(:,1) +  halp(:,3) .* hbeta(:,2)+halp(:,4).*hbeta(:,3);
 res = y-fit;
 p = 4*(m+k1)+3*(m+k2);
 sig = sum(res.^2)/(T-p);
end

