function [hbeta,alp0,fit,res,sig, coeff, c, up0]= RM( sampleX,u,y,m,k1,k2,delta )
%construct spline estimation of additive model with selected three
%covariates
T=length(y);
%B-spline matrix
Bt = rspline(u,u,m,k1);
B1 = rspline(sampleX(:,1),sampleX(:,1),m,k2);
B2 = rspline(sampleX(:,2),sampleX(:,2),m,k2);
B3 = rspline(sampleX(:,3),sampleX(:,3),m,k2);
B = [Bt B1 B2 B3];
coeff =  pinv(B'*B+delta*eye(size(B,2)))*B'*y;
%spline estimation of additive function
BX= blkdiag(Bt, B1, B2, B3);
est = BX*coeff;
est = reshape(est, T,4);
c = mean(est(:,2:end));
hbeta = est(:, 2 : end)-repmat(c,T,1);
alp0 = est(:, 1);
up0 = c(1)+c(2)+c(3);
fit = up0 + alp0+hbeta(:,1)+hbeta(:,2)+hbeta(:,3);
res = y-fit;
sig =sum(res.^2)/(T - 3*m-3*k2-m-k1); 
end


