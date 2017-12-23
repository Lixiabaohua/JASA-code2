function [res, fit, coeff1, coeff2, a, c] =RM( Inib, k1, k2, m, sampleX, u, y, delta )
% resduced model estimation
T=length(y);
% Step II estimation
 Bt= rspline(u, u ,m, k1); 
 D1=(kron(Inib(:,1), ones(1,size(Bt,2)))).*Bt;
 D=[Bt D1 Inib(:,2)];
 coeff1  = pinv(D' * D+delta * eye(size(D,2)))*D' * y;
 %estimation fot varying-coefficient function
 halp0 = Bt *  coeff1(1: (k1 + m));
 halp1 = Bt *  coeff1((k1 + m + 1) : 2 * (k1 + m));
 halp2 = coeff1( 2 * (k1 + m) + 1);
 %standardized
 a = mean(halp1) ;
 halp1 = halp1/a; 
 % Step III B-spline  to estimate additive functions
 esty = y - halp0;
 Bx1 = rspline(sampleX(:,1), sampleX(:,1), m, k2);   
 Bx2 = rspline(sampleX(:,2), sampleX(:,2), m, k2);   
 S1 =  (kron(halp1, ones(1,  size(Bx1,2)))) .* Bx1;
 S   = [S1 halp2 * Bx2];
 coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S' * esty;
 %estimation for additive functions
 BX= blkdiag(Bx1,Bx2);
 hbeta = BX *coeff2;
 hbeta = reshape(hbeta,T,2);
 %centering
 c =mean(hbeta);
 hbeta = hbeta-repmat(c,T,1);
 %updated estimation for alpha0
 up0 = halp0 +  halp1 * c(1) + halp2 * c(2);
 fit = up0 + halp1.* hbeta(:,1) + halp2 * hbeta(:,2);
 res = y-fit;
end





