function [ estMat,coeff1,coeff2,mynorm,c] = SpestP(y,b1,b2,Bt,Bx1,Bx2,I1,delta )
%  I1: segment length in step I estimation
% Bx1,Bx2,Bt,b1,b2:B-spline basis matrix
% delta: ridge regression parameters
% coeff1: coefficients for varying-coefficient functions
% coeff2: coefficients for additive functions

%sample size
T=length(y);   
N = round(T/I1)-1;
%estimation in Step I
B = [ones(T,1) b1 b2 ];
s1 = zeros(size(b1,2),1);  s2 = zeros(size(b2,2),1);  
for l = 1 : N
      low = (l - 1 ) *I1 +1 ;
      up = I1* l ;
      bx = B(low : up, : ) ;
      by = y( low : up) ;
      tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
      s1 = s1 + tem3(2:size(b1,2)+1);
      s2 = s2 + tem3(size(b1,2)+2:2*size(b1,2)+1);
end  
low = N*I1 + 1;
up =T;
bx = B(low : up, : ) ;
by = y( low : up) ;
tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
s1 = s1 + tem3(2:size(b1,2)+1);
s2 = s2 + tem3(size(b1,2)+2:2*size(b1,2)+1);
%initial estimate for additive function
hbeta1 =  b1 * s1;
hbeta2 =  b2 * s2;
%centering
Inibeta1 = hbeta1 -mean(hbeta1);
Inibeta2 = hbeta2 -mean(hbeta2);
%Estimation in  Step II 
D1=(kron(Inibeta1, ones(1,size(Bt,2)))).*Bt;
D2=(kron(Inibeta2, ones(1,size(Bt,2)))).*Bt;
D=[Bt D1 D2];
coeff1  = pinv(D'*D+delta*eye(size(D,2)))*D'*y;
%estimation fot varying-coefficient function
B = blkdiag(Bt,Bt, Bt);
alpvec = B*coeff1;
ealp = reshape(alpvec,T,3);
%standard 
mynorm = mean(ealp(:,2:3)) ;
ealp = ealp./repmat([1 mynorm],T,1);
 % Estimation in Step III
pseuy =y - ealp(:,1) ;
S1 =  (kron(ealp(:,2), ones(1,size(Bx1,2)))).*Bx1;
S2 =  (kron(ealp(:,3), ones(1,size(Bx2,2)))).*Bx2;
S   = [S1 S2 ];
coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S'*pseuy;
BX = blkdiag(Bx1,Bx2);
betavec = BX*coeff2;
ebeta = reshape(betavec,T,2);
c = mean(ebeta);
ealp(:,1) = ealp(:,1) + ealp(:,2) *c(1)+ealp(:,3)*c(2);
estMat = [ealp  ebeta];
end

