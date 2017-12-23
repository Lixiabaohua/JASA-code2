function  [err, yA]= RM_test( T, X, u,y,n,m,k1,k2,delta)
%return one-step  prediction and prediction erro
Bt = rspline(u,u,m,k1);
Bx1 = rspline(X(:,1),X(:,1), m,k2);   
Bx2 = rspline(X(:,2),X(:,2), m,k2); 
Bx3 = rspline(X(:,3),X(:,3), m,k2); 
yA =  zeros(n,1);
for i =1:n
Acoeff = pre_RM(y(1:T-n+i-1),Bt(1:T-n+i-1,:),...
                         Bx1(1:T-n+i-1,:),Bx2(1:T-n+i-1,:), Bx3(1:T-n+i-1,:),delta ) ;
BX = blkdiag(Bt(T-n+i,:),Bx1(T-n+i,:),Bx2(T-n+i,:),Bx3(T-n+i,:));
est = BX * Acoeff;
%beta = est(2:end)-c';
%alp0 = est(1)+c(1)+c(2)+c(3);
yA(i) = est(1)+ est(2) + est(3) + est(4)  ;
end
err = y(T-n+1:end) -yA;
end

