function [ alp0, beta] = pred_RMest(px1, px2, px3, pt, sampleX, u, coeff, c, up0, m, kC, kA )
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
MatX3 = rspline(sampleX(:, 3), px3, m, kA);
alp0 = MatT * coeff(1 : size(MatT,2)) + up0; 
%% prediction for additive function
beta1 = MatX1 * coeff((size(MatT, 2)+1) : (size(MatT, 2) + size(MatX1, 2)));
beta2 = MatX2 * coeff((size(MatT, 2) + size(MatX1, 2) + 1) : (size(MatT, 2) + 2 * size(MatX1, 2) ));
beta3 = MatX3 * coeff((size(MatT, 2) + 2*size(MatX1, 2) + 1) : (size(MatT, 2) + 3 * size(MatX1, 2) ));
beta = [beta1 beta2 beta3] -repmat(c, size(MatX1, 1), 1);
end




