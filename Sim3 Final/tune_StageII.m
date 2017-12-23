function [const,Pbeta,Pnbeta,mybic,Pnline]= tune_StageII(Palp,k,m,mu,sampleX,y,bound,p,delta)
%choose tuning parameter in stage II and get penalized estimation of additive functions
% Palp: penalized estimation of varying-coefficient function
% mu: tuning parameter in Stage II 
% sampleX,y: observations
% k: knots for additive function
% m: order of B-spline
% p: dimension of covariates
% delta : ridge regression parameter
% bound: convergence bound

T=length(y);
%additive function B-spline
[Bx1,DBx1] = rspline2(sampleX(:,1), sampleX(:,1),m,k); 
[Bx2,DBx2] = rspline2(sampleX(:,2), sampleX(:,2),m,k); 
[Bx3,DBx3] = rspline2(sampleX(:,3), sampleX(:,3),m,k); 
[Bx4,DBx4] = rspline2(sampleX(:,4), sampleX(:,4),m,k);
% to find 
[~,D] = rspline(sampleX(:,4), sampleX(:,4),m,k); 
%weight matrix in stage II model identification
V1 = DBx1'*DBx1/T;  V2 = DBx2'*DBx2/T;
V3 = DBx3'*DBx3/T;  V4 = DBx4'*DBx4/T;
Z = D'*D/T;      
%substuting into penalized varying-coefficient estimation
Py = y - Palp(:,1);
SP1 = (kron(Palp(:,2), ones(1,size(Bx1,2)))) .* Bx1;
SP2 = (kron(Palp(:,3), ones(1,size(Bx2,2)))) .* Bx2;
SP3 = (kron(Palp(:,4), ones(1,size(Bx3,2)))) .* Bx3;
SP4 = (kron(Palp(:,5), ones(1,size(Bx4,2)))) .* Bx4;
SP=[SP1 SP2 SP3 SP4];
Pcoeff2 = pinv((SP)'*SP+delta * eye(size(SP,2)))*(SP)'*Py; 
diff2 = 1; coeff2 = Pcoeff2;
Omega2 = zeros(size(SP,2));
while (diff2>bound)
 coefind1 = coeff2(1:size(Bx1,2));
 coefind2 = coeff2(size(Bx1,2)+1:2*size(Bx1,2));
 coefind3 = coeff2(2*size(Bx1,2)+1:3*size(Bx1,2));
 coefind4 = coeff2(3*size(Bx1,2)+1:4*size(Bx1,2));
 nbeta1 =  sqrt(abs(coefind1'*V1*coefind1));
 nbeta2 =  sqrt(abs(coefind2'*V2*coefind2));
 nbeta3 =  sqrt(abs(coefind3'*V3*coefind3));
 nbeta4 =  sqrt(abs(coefind4'*V4*coefind4));
 constbeta1 = deriveSCAD(k^(0)*nbeta1,mu)/nbeta1; 
 constbeta2 = deriveSCAD(k^(0)*nbeta2,mu)/nbeta2; 
 constbeta3=  deriveSCAD(k^(0)*nbeta3,mu)/nbeta3;
 constbeta4=  deriveSCAD(k^(0)*nbeta4,mu)/nbeta4;
 Omega2(1:size(Bx1,2),1:size(Bx1,2))= constbeta1*V1;
 Omega2(size(Bx1,2)+1:2*size(Bx1,2),size(Bx1,2)+1:2*size(Bx1,2))= constbeta2*V2;
 Omega2(2*size(Bx1,2)+1: 3*size(Bx1,2),2*size(Bx1,2)+1: 3*size(Bx1,2)) = constbeta3*V3;
 Omega2(3*size(Bx1,2)+1: 4*size(Bx1,2),3*size(Bx1,2)+1: 4*size(Bx1,2)) = constbeta4*V4;
 tem = pinv((SP)'*SP+T*Omega2+delta*eye(size(SP,2)))*(SP)'*Py;
 diff2 = norm(tem-coeff2);
 coeff2 = tem; 
end
 %penalized estimation of additive function \beta_{k}
coefind1 = coeff2(1:size(Bx1,2));
coefind2 = coeff2(size(Bx1,2)+1:2*size(Bx1,2));
coefind3 = coeff2(2*size(Bx1,2)+1:3*size(Bx1,2));
coefind4 = coeff2(3*size(Bx1,2)+1:4*size(Bx1,2));
Pbeta1 = Bx1 * coefind1;
Pbeta2 = Bx2 * coefind2;
Pbeta3 = Bx3 * coefind3;
Pbeta4 = Bx4 * coefind4;
Pbeta = [Pbeta1 Pbeta2 Pbeta3 Pbeta4];
c = mean(Pbeta);
Pbeta = Pbeta - repmat(c,T,1);
const = Palp(:,1) + Palp(:,2)*c(1)+Palp(:,3)*c(2)+Palp(:,4)*c(3) +Palp(:,5)*c(4);
Pnbeta1 = sqrt(abs(coefind1'*V1*coefind1));
Pnbeta2 = sqrt(abs(coefind2'*V2*coefind2));
Pnbeta3 = sqrt(abs(coefind3'*V3*coefind3));
Pnbeta4 = sqrt(abs(coefind4'*V4*coefind4));
Pnline= sqrt(abs(coefind4'*Z*coefind4));
Pnbeta = [Pnbeta1 Pnbeta2 Pnbeta3 Pnbeta4 ];
Ind = find(Pnbeta<bound);   %varying-coefficient term index
varyterm = length(Ind);
fit =const + Palp(:,2).*Pbeta(:,1)+Palp(:,3).*Pbeta(:,2) + Palp(:,4).*Pbeta(:,3)+...
        Palp(:,5).*Pbeta(:,4) ;
rss2 = sum((y- fit ).^2);    
mybic = log(rss2)+(varyterm + (p-varyterm)*(m+k))* log(T)/T;
end

