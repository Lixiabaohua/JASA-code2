function [fit, res] = BivEst( sampleX, u, y, m, kC, kA,delta)
% construct estimations of bivariate additive component functions 
% sampleX is matrix of covariates
bx1 = rspline(sampleX(:,1), sampleX(:,1), m, kA);   
bx2 = rspline(sampleX(:,2), sampleX(:,2), m, kA); 
bt= rspline(u, u, m, kC); 
% tensor product of B-spline 
tensb1 = kron(bx1,ones(1,size(bt,2))) .* repmat(bt,1,size(bx1,2));
tensb2 = kron(bx2,ones(1,size(bt,2))) .* repmat(bt,1,size(bx2,2));
xMat = [bt tensb1 tensb2] ;
coeff = pinv(xMat'*xMat + delta*eye(size(xMat,2)))*xMat'*y;
fit = bt * coeff(1 : size(bt,2)) + tensb1 * coeff((1 + size(bt,2)) : (size(tensb1,2) + size(bt,2)))...
           +  tensb2*coeff((size(tensb1,2) + size(bt,2)+1): (2*size(tensb1,2) + size(bt,2)));
res = y - fit;
end

