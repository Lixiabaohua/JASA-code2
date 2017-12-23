function [halp,hbeta,res] =Spest( I1,kC, kA, m1, m, sampleX, u, y,delta)
%Step II and Step III estimation (using same order of B-spline
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline in Step II and III estimation
% m1: order of B-spline in Step I estimation
% I1: segment length in step I estimation
% sampleX,u,y are observations
%delta: ridge regression parameter

%initial estimation in Step I
T=length(y);
b1 = rspline(sampleX(:,1),sampleX(:,1), m1,kA);   
b2 = rspline(sampleX(:,2),sampleX(:,2), m1,kA);
b3 = rspline(sampleX(:,3),sampleX(:,3), m1,kA);
b4 = rspline(sampleX(:,4),sampleX(:,4), m1,kA);   
%design matrix in Step I estiamtion
B = [ones(T,1) b1 b2 b3 b4];
%Step I estimation: segment 
s1 = zeros(size(b1,2),1);  s2 = zeros(size(b2,2),1);  
s3 = zeros(size(b3,2),1);  s4 = zeros(size(b4,2),1);  
N = T/I1;  
for l = 1 : N
      low = (l - 1 ) * I1 +1 ;
      up = I1 * l ;
      bx = B(low : up, : ) ;
      by = y( low : up) ;
      tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
      s1 = s1 + tem3(2:size(b1,2)+1);
      s2 = s2 + tem3(size(b1,2)+2:2*size(b1,2)+1);
      s3 = s3 + tem3(2*size(b1,2)+2:3*size(b1,2)+1);
      s4 = s4 + tem3(3*size(b1,2)+2:4*size(b1,2)+1);
end       
%initial estimate for additive function
hbeta1 =  b1 * s1;
hbeta2 =  b2 * s2;
hbeta3 =  b3 * s3;
hbeta4 =  b4 * s4;
%centering
hbeta = [hbeta1 hbeta2 hbeta3 hbeta4];
Ini =  hbeta-repmat(mean(hbeta),T,1);
% estion of Step II (varying-coefficient functions)
Bt= rspline(u, u ,m, kC); 
D1 = (kron(Ini(:,1), ones(1,size(Bt,2)))).*Bt;
D2 = (kron(Ini(:,2), ones(1,size(Bt,2)))).*Bt;
D3 = (kron(Ini(:,3), ones(1,size(Bt,2)))).*Bt;
D4 = (kron(Ini(:,4), ones(1,size(Bt,2)))).*Bt;
D=[Bt D1 D2 D3 D4];
coeff1    = pinv(D'*D+delta*eye(size(D,2)))*D'*y;
B= blkdiag(Bt,Bt,Bt,Bt,Bt);
halp = B*  coeff1;
halp = reshape(halp,T,5);
%standardized
mynorm = mean(halp(:,2:end));
halp = halp./repmat([1 mynorm],T,1);
% estimation of Step III (additive functions)
Bx1 = rspline(sampleX(:,1),sampleX(:,1), m,kA); 
Bx2 = rspline(sampleX(:,2),sampleX(:,2), m,kA); 
Bx3 = rspline(sampleX(:,3),sampleX(:,3), m,kA); 
Bx4 = rspline(sampleX(:,4),sampleX(:,4), m,kA); 
pseuy =y - halp(:,1) ;
S1 =  (kron(halp(:,2), ones(1,size(Bx1,2)))).*Bx1;
S2 =  (kron(halp(:,3), ones(1,size(Bx2,2)))).*Bx2;
S3 =  (kron(halp(:,4), ones(1,size(Bx3,2)))).*Bx3;
S4 =  (kron(halp(:,5), ones(1,size(Bx4,2)))).*Bx4;
S   = [S1 S2 S3 S4];
coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S'*pseuy;
%estimation for additive functions
BX= blkdiag(Bx1,Bx2,Bx3, Bx4);
hbeta = BX *coeff2;
hbeta = reshape(hbeta,T,4);
%centering for additive function
c  = mean(hbeta);
hbeta = hbeta - repmat(c,T,1);
halp(:,1) = halp(:,1) +halp(:,2)*c(1) +halp(:,3)*c(2) +halp(:,4)*c(3)+...
                halp(:,5)*c(4);
fit = halp(:,1) + halp(:,2).* hbeta(:,1) +  halp(:,3) .* hbeta(:,2)+...
        halp(:,4).*hbeta(:,3) + halp(:,5).*hbeta(:,4);
res = y-fit;
end

