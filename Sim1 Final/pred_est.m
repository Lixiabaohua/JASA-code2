function [ halp,hbeta ] = pred_est(px1,px2,pt,sampleX,u,coeffC,coeffA,stpara,m,kC,kA )
% predict at given grid
% px1,px2,pt: given grid points
% prediction for varying-coefficient function
% sampleX,u: sample
% coeffC,coeffA: spline coefficients estimationg varying-coefficient and
% additive function
% stpara: standard parameter giving norm of varying-coefficient function
% and center of additive function
% order of B-spline function
% kC, kA: interior knots numner
  MatT = rspline(u, pt,m,kC);
MatX1 = rspline(sampleX(:,1),px1,m,kA);
MatX2 = rspline(sampleX(:,2),px2,m,kA);
B= blkdiag(MatT,MatT,MatT);
halp = B*coeffC;
halp = reshape(halp,length(pt),3);
%standardized
halp = halp./repmat([1 stpara(1:2)],length(pt),1);
%estimationprediction for additive function
BX = blkdiag(MatX1,MatX2);
hbeta = BX*coeffA;
hbeta = reshape(hbeta,length(px1),2);
hbeta = hbeta - repmat(stpara(3:4),length(px1),1);
end

