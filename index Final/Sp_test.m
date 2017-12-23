function [err, ypre] = Sp_test(T, X, u, y, n, kC, kA, m1, m, delta)
%compute the estimation on traing set and prediction on testing test
% n: length of prediction
% T: sample length 
%Step II
Bt= rspline(u, u ,m, kC); 
%Step III
Bx1 = rspline(X(:,1),X(:,1), m,kA);   
Bx2 = rspline(X(:,2),X(:,2), m,kA); 
ypre = zeros(n,1);
for i =1:n
  Ini= StepIpre( kA, m1, X(1:T-n+i-1,:),  y(1:T-n+i-1), delta ) ;
  [ ~, ~, coefC, coefA, a] = RM( Ini, kC, kA, m, X(1:T-n+i-1,:), u(1:T-n+i-1), y(1:T-n+i-1), delta ) ;
  alpvec1 = Bt(T-n+i,:) * coefC(1: (kC+m));
  alpvec2 = Bt(T-n+i,:) * coefC((kC + m + 1) : (2 * (kC + m)));
  alpvec3 =  coefC(2 * (kC + m)  + 1);
  alpvec2 = alpvec2/a;
 BX = blkdiag(Bx1(T-n+i,:),Bx2(T-n+i,:));
 betavec = BX * coefA;
 ypre(i) = alpvec1+ alpvec2 *betavec(1)+...
            alpvec3 * betavec(2);
end
err = y(T-n+1 : end) - ypre;
end

