function [ alp,beta] = pred_RMest(px1, px2, pt, sampleX, u, coeffC, coeffA, a, c, m, kC, kA )
%predict of component functions at given grid
% px1, px2 and pt are prediction points
% sampleX and u are sample observations
% coeffC,coeffA are coefficients estimating component functions
% a and c are norm and mean of component functions
% kC and kA are optimal knots
% m: order of b-spline
%% prediction for varying-coefficient function
  MatT = rspline(u, pt, m, kC);
MatX1 = rspline(sampleX(:, 1), px1, m, kA);
MatX2 = rspline(sampleX(:, 2), px2, m, kA);
halp0 = MatT * coeffC(1 : (m + kC)); 
halp1 = MatT * coeffC((m + kC +1) : 2 * (m + kC)) ;
alp1 = halp1/a;
%% prediction for additive function
hbeta1 = MatX1 * coeffA(1 : (m+kA));
beta1 = hbeta1 - c(1);
hbeta2 = MatX2*coeffA((m + kA + 1) :  2 * (m + kA));
beta2 = hbeta2 - c(2);
alp = [halp0 alp1];
beta = [beta1 beta2];
end


