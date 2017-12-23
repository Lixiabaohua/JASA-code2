function [Acoeff,c] = pre_RM(y,Bt,Bx1,Bx2,Bx3,delta)
%additive model estimation for one-step ahead forecasting
T = length(y);
Bx = [Bt Bx1 Bx2 Bx3];
Acoeff =  pinv(Bx'*Bx+delta*eye(size(Bx,2)))*Bx'*y;
BX = blkdiag(Bt,Bx1,Bx2,Bx3);
betavec = BX*Acoeff;
ebeta = reshape(betavec,T,4);
c = mean(ebeta(:,2:end));
end

