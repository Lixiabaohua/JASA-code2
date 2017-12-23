function rory= Sp_test(T, Dat, repy, n, kC,kA,m1,m,gnum,delta )
%compute one-step ahead forecasting based on VCAM
% n: length of prediction
% gnum: segment length in Step I
% T: sample length
% kC: knots for varying-coefficient function
% kA:  knots for additive function
% m1, m: order of B-spline 
% delta: ridge regression parameter
% repy: Monte Carlo replication of response
% Dat: data set generated
%Step I 
ory = Dat(:,1)+Dat(:,2).*Dat(:,4) + Dat(:,3).*Dat(:,5); %mean response
b1 = rspline(Dat(:,6),Dat(:,6), m1,kA);   
b2 = rspline(Dat(:,7),Dat(:,7), m1,kA);  
%Step II
Bt= rspline(Dat(:,8), Dat(:,8) ,m, kC); 
%Step III
Bx1 = rspline(Dat(:,6),Dat(:,6), m,kA);   
Bx2 = rspline(Dat(:,7),Dat(:,7), m,kA); 
ypre = zeros(n,1);
for i = 1:n
[ ~,coeffC,coeffA,stdpara,c] = SpestP(repy(1:T-n+i-1),b1(1:T-n+i-1,:),...
    b2(1:T-n+i-1,:),Bt(1:T-n+i-1,:),Bx1(1:T-n+i-1,:),Bx2(1:T-n+i-1,:),gnum,delta ) ;
B = blkdiag(Bt(T-n+i,:),Bt(T-n+i,:),Bt(T-n+i,:));
alpvec  = B *coeffC;
alp = alpvec./[1 stdpara]' ;
BX = blkdiag(Bx1(T-n+i,:),Bx2(T-n+i,:));
beta = BX * coeffA;
beta = beta - c';
alp(1) = alp(1)+alp(2)*c(1)+alp(3)*c(2);
ypre(i) = alp(1) + alp(2).*beta(1)+alp(3).*beta(2);
end
rory = sqrt(mean((ypre-ory(T-n+1:end)).^2));
end

