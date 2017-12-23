function [hbeta,const,fit,res,sig]= add_est( sampleX,y,m,k,delta )
%construct spline estimation of additive model with selected three
%covariates
T=length(y);
%B-spline matrix
B1 = rspline(sampleX(:,1),sampleX(:,1),m,k);
B2 = rspline(sampleX(:,2),sampleX(:,2),m,k);
B = [ones(T,1) B1 B2];
coeff =  pinv(B'*B+delta*eye(size(B,2)))*B'*y;
%spline estimation of additive function
BX= blkdiag(B1,B2);
hbeta = BX*coeff(2:end);
hbeta = reshape(hbeta, T,2);
c = mean(hbeta);
hbeta = hbeta-repmat(c,T,1);
const =coeff(1)+ c(1)+c(2);
fit = const+hbeta(:,1)+hbeta(:,2);
res = y-fit;
sig =sum(res.^2)/(T - 2*m-2*k-1); 
end


