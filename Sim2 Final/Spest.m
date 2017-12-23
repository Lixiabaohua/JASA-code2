function [fit,res] =Spest( I1,kC, kA, m1, m, sampleX, u, y,delta)
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline in Step II and III estimation
% m1: order of B-spline in Step I estimation
% I1: segment length in step I estimation
% sampleX, u, y:sample
% delta: redge regression parameter
%initial estimation in Step I
T=length(y);
b1 = rspline(sampleX(:,1), sampleX(:,1), m1,kA);   
b2 = rspline(sampleX(:,2), sampleX(:,2), m1,kA); 
%design matrix in Step I estiamtion
B = [ones(T,1) b1 b2 ];
%Step I estimation: segment 
s1 = zeros(size(b1,2),1);  s2 = zeros(size(b2,2),1);  
N =T/I1;
for l = 1 : N
      low = (l - 1 ) * I1 +1 ;
      up = I1 * l ;
      bx = B(low : up, : ) ;
      by = y( low : up) ;
      tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
      s1 = s1 + tem3(2:(size(b1,2)+1));
      s2 = s2 + tem3((size(b1,2)+2):(2*size(b1,2)+1));
end       
%initial estimate for additive function
hbeta1 =  b1 * s1;
hbeta2 =  b2 * s2;
%centering
Inibeta1 = hbeta1 - mean(hbeta1);
Inibeta2 = hbeta2 - mean(hbeta2);
% estion of Step II (varying-coefficient functions)
Bt= rspline(u, u ,m, kC); 
D1 = (kron(Inibeta1, ones(1,size(Bt,2)))).*Bt;
D2=(kron(Inibeta2, ones(1,size(Bt,2)))).*Bt;
D=[Bt D1 D2];
coeff1    = pinv(D'*D+delta*eye(size(D,2)))*D'*y;
B= blkdiag(Bt, Bt, Bt);
halp = B*  coeff1;
halp = reshape(halp,T,3);
%standardized
mynorm = mean(halp(:,2:3));
halp = halp./repmat([1 mynorm],T,1);
% estimation of Step III (additive functions)
Bx1 = rspline(sampleX(:,1),sampleX(:,1), m,kA);   
Bx2 = rspline(sampleX(:,2),sampleX(:,2), m,kA);
pseuy =y - halp(:,1) ;
S1 =  (kron(halp(:,2), ones(1,size(Bx1,2)))).*Bx1;
S2 =  (kron(halp(:,3), ones(1,size(Bx2,2)))).*Bx2;
S   = [S1 S2 ];
coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S'*pseuy;
%estimation for additive functions
BX= blkdiag(Bx1,Bx2);
hbeta = BX *coeff2;
hbeta = reshape(hbeta,T,2);
c = mean(hbeta);
hbeta = hbeta - repmat(c,T,1);
halp(:,1) = halp(:,1) + halp(:,2)*c(1)+halp(:,3)*c(2);
fit = halp(:,1) + halp(:,2).* hbeta(:,1) +  halp(:,3) .* hbeta(:,2);
res = y-fit;
%stdpara = [mynorm c];
end

