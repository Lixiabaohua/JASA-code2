 function [Palp,Pnalp,mybic]= tune_StageI( funbeta,k1,m,lam,u,y,bound,p,delta)
% Stage I penalized estimation to identify constant coefficient function
% get penalized function estimation and its norm of the first derivative of varying-coefficient function
% k1 for varying-coefficient function
% m: order of B-spline
% lam:tuning parameter for varying-coefficient functions
% u,y :observations
% funbeta: spline estimation of additive functions
% p: dimension of covariates
% delta : ridge regression parameter
% bound: convergence bound

[Bt,B1t] = rspline(u,u,m,k1); %Step II
T=length(y);
 %penalizing matrix
W = (B1t)'*B1t/T;          
%substituing spline estimation of beta function 
DP1=(kron(funbeta(:,1), ones(1,m+k1))).*Bt;  
DP2=(kron(funbeta(:,2), ones(1,m+k1))).*Bt;   
DP3=(kron(funbeta(:,3), ones(1,m+k1))).*Bt;
DP4=(kron(funbeta(:,4), ones(1,m+k1))).*Bt;   
DP=[Bt DP1 DP2 DP3 DP4];
Pcoeff1 = pinv((DP)'*DP+delta*eye(size(DP,2)))*(DP)'*y; 
%initial value of while looping
diff1 = 1;  coeff1 =Pcoeff1;
 while (diff1>bound)
 coefind1 = coeff1(size(Bt,2)    + 1 : 2*size(Bt,2));
 coefind2 = coeff1(2*size(Bt,2) + 1 : 3*size(Bt,2));
 coefind3 = coeff1(3*size(Bt,2) + 1 : 4*size(Bt,2));
 coefind4 = coeff1(4*size(Bt,2) + 1 : 5*size(Bt,2));
 temp1 = Bt*coefind1;
 temp2 = Bt*coefind2;
 temp3 = Bt*coefind3;
 temp4 = Bt*coefind4;
 %norm of first derivative of estimated varying-coefficient (\alpha_{k}) function
 nalpha1 = sqrt(abs(coefind1'*W*coefind1))/mean(temp1);
 nalpha2 = sqrt(abs(coefind2'*W*coefind2))/mean(temp2);
 nalpha3 = sqrt(abs(coefind3'*W*coefind3))/mean(temp3);
 nalpha4 = sqrt(abs(coefind4'*W*coefind4))/mean(temp4);
 constalp1 = deriveSCAD(k1^(0)*nalpha1,lam)/nalpha1; 
 constalp2 = deriveSCAD(k1^(0)*nalpha2,lam)/nalpha2; 
 constalp3=  deriveSCAD(k1^(0)*nalpha3,lam)/nalpha3;
 constalp4=  deriveSCAD(k1^(0)*nalpha4,lam)/nalpha4; 
 Omega1 = kron(diag([0,constalp1,constalp2,constalp3,constalp4]),W);
 tem = pinv((DP)'*DP+T*Omega1+delta*eye(size(DP,2)))*(DP)'*y;
 diff1 = norm(tem-coeff1);
 coeff1 = tem; 
 end
 %penalized estimation of varying-coefficient function
 coefind0 = coeff1(1:size(Bt,2));
 coefind1 = coeff1(size(Bt,2)+1:2*size(Bt,2));
 coefind2 = coeff1(2*size(Bt,2)+1:3*size(Bt,2));
 coefind3 = coeff1(3*size(Bt,2)+1:4*size(Bt,2));
 coefind4 = coeff1(4*size(Bt,2)+1:5*size(Bt,2));
 iterPalp0=Bt*coefind0;
 iterPalp1=Bt*coefind1;
 iterPalp2=Bt*coefind2;
 iterPalp3=Bt*coefind3;
 iterPalp4=Bt*coefind4;
 a1 = mean(iterPalp1);
 a2 = mean(iterPalp2);
 a3 = mean(iterPalp3);
 a4 = mean(iterPalp4);
 %standardized
 iterPalp1 =  iterPalp1/a1;
 iterPalp2 =  iterPalp2/a2;
 iterPalp3 =  iterPalp3/a3;
 iterPalp4 =  iterPalp4/a4;
 Pnalpha1 = sqrt(abs(coefind1'*W*coefind1))/a1;
 Pnalpha2 = sqrt(abs(coefind2'*W*coefind2))/a2;
 Pnalpha3 = sqrt(abs(coefind3'*W*coefind3))/a3;
 Pnalpha4 = sqrt(abs(coefind4'*W*coefind4))/a4;
 Pnalp = [ Pnalpha1 Pnalpha2 Pnalpha3 Pnalpha4 ];
 invInd = find(Pnalp<bound);     %invariant coefficient function index
 addterm = length(invInd);         %number of  invariant coefficient function 
 fit = iterPalp0 + iterPalp1.*funbeta(:,1) +iterPalp2.*funbeta(:,2) + iterPalp3.*funbeta(:,3)+...
         iterPalp4.*funbeta(:,4);
 rss1 = sum((y-fit ).^2);
 mybic= log(rss1)+(addterm +(p+1-addterm)*(m+k1))* log(T)/T;
 Palp = [iterPalp0 iterPalp1 iterPalp2 iterPalp3 iterPalp4];
end

