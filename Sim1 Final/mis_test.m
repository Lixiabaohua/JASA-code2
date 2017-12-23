function [ryVn,ryAn]= mis_test( T, Dat, y, n, kC,kA,m,delta)
%compute one-step ahead forecasting under misspecified VCM and AM
% n: length of prediction
% T: sample length
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m: order of B-spline 
% delta: ridge regression parameter
% y: Monte Carlo replication of response
% Dat: data set generated

ory = Dat(:,1)+Dat(:,2).*Dat(:,4) + Dat(:,3).*Dat(:,5); %response without error
Bt= rspline(Dat(:,8), Dat(:,8) ,m, kC); 
Bx1 = rspline(Dat(:,6),Dat(:,6), m,kA);   
Bx2 = rspline(Dat(:,7),Dat(:,7), m,kA); 
yV = zeros(n,1); yA = zeros(n,1);
for i = 1:n
[~,~,Vcoeff,Acoeff,conpara ] = mis_Pest(Dat(1:T-n+i-1,6:7),y(1:T-n+i-1),Bt(1:T-n+i-1,:),...
                         Bx1(1:T-n+i-1,:),Bx2(1:T-n+i-1,:),delta ) ;
B =blkdiag(Bt(T-n+i,:),Bt(T-n+i,:),Bt(T-n+i:end,:));
alp  = B *Vcoeff;
BX = blkdiag(Bx1(T-n+i,:),Bx2(T-n+i,:));
betavec = BX * Acoeff(2:end);
beta= [conpara(3) betavec' - conpara(1:2)];
yV(i) = alp(1)+ alp(2).*Dat(T-n+i,6)+alp(3).*Dat(T-n+i,7);
yA(i) = beta(1) + beta(2) +beta(3) ;
end
ryVn = sqrt(mean((yV-ory(T-n+1:end)).^2));
ryAn = sqrt(mean((yA-ory(T-n+1:end)).^2));
end

