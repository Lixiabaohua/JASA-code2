    
function  [Bt, B1t] = rspline(t,t0,m,k)
% using sample t to generate b-spline basis matrix with (m-1)-th order as well as its 
% deravative matrix of given points t0 with k uniform internal knots on t
% % Attention! it must be gauranteed that there is no replicated
% samples !!!!!!!
[r1] = size(t0,1);
[tt,order] = sort(t0);
tt1 = zeros(k,1);
for i=1:k
    tt1(i,1)=prctile(sort(t),i*100/(k+1));
end
b = max(t) + 10^(-10);
a = min(t) - 10^(-10);
tt2=[a*ones(1,m),tt1',b*ones(1,m)]';
ttt = kron(tt,ones(2,1)); %每一行重复两次
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Amatrix = spcol(tt2,m,ttt);
A1 = Amatrix(1:2:2*r1-1,:);
A2 = Amatrix(2:2:2*r1,:);
Bt = zeros(r1,k+m);
B1t = zeros(r1,k+m);
for i =1:r1
    Bt(order(i),:) = A1(i,:);
    B1t(order(i),:) = A2(i,:);
end
end



