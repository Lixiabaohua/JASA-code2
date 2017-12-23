function [ihalp,ihbeta,res,coeff1,coeff2,stdpara] =IterSpest(sampleX,u,y,hbeta,kA,kC,m,delta)
% step I: iterative varying-coefficient component function
% step II: iterative additive component function
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline in Step II and III estimation
% sampleX, u, y:sample
% delta: redge regression parameter
%initial estimation in Step I
T=length(y);
 
% step I: iterative varying-coefficient component function
Bt= rspline(u, u ,m, kC); 
D1 = (kron(hbeta(:,1), ones(1,size(Bt,2)))).*Bt;
D2=(kron(hbeta(:,2), ones(1,size(Bt,2)))).*Bt;
D=[Bt D1 D2];
coeff1  = pinv(D'*D+delta*eye(size(D,2)))*D'*y;
B= blkdiag(Bt,Bt,Bt);
ihalp = B*  coeff1;
ihalp = reshape(ihalp,T,3);
%standardized
mynorm = mean(ihalp(:,2:3));
ihalp = ihalp./repmat([1 mynorm],T,1);

% step II: iterative additive component function
Bx1 = rspline(sampleX(:,1),sampleX(:,1), m,kA);   
Bx2 = rspline(sampleX(:,2),sampleX(:,2), m,kA);
pseuy =y - ihalp(:,1) ;
S1 =  (kron(ihalp(:,2), ones(1,size(Bx1,2)))).*Bx1;
S2 =  (kron(ihalp(:,3), ones(1,size(Bx2,2)))).*Bx2;
S   = [S1 S2 ];
coeff2 =  pinv(S'*S+delta*eye(size(S,2)))*S'*pseuy;

%estimation for additive functions
BX= blkdiag(Bx1,Bx2);
ihbeta = BX *coeff2;
ihbeta = reshape(ihbeta,T,2);
c = mean(ihbeta);
%centralized
ihbeta = ihbeta - repmat(c,T,1);

ihalp(:,1) = ihalp(:,1) + ihalp(:,2)*c(1)+ihalp(:,3)*c(2);

fit = ihalp(:,1) + ihalp(:,2).* ihbeta(:,1) +  ihalp(:,3) .* hbeta(:,2);
res = y-fit;
stdpara = [mynorm c];
end


