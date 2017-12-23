function [ Msig2,ind ] = estStd(Dat,Q,kA,kC,m,I,delta)
%check the performance of standard deviation(based on 1000 Monte Carlo
%replications)
% I,kA,kC are selected segment length and knot number
% output median estimates of variance of rendom error

% sample in estimating
T=size(Dat,1);
t = Dat(:,8);
X=Dat(:,6:7);

%B-spline approximating 
Bt= rspline(t, t ,m, kC); 

My=MontY(Dat,Q);

% true value of  standard deviation and variance of random error 
truesig=0.4*(0.7+(2-t)./(2+t));
truesig2=truesig.^2;

%Mse restore mse in each replication; sig2 restore the estimated variance
%of random error
Mse=zeros(1,Q);  sig2=zeros(T,Q);

 for i = 1: Q
    %three-step spline estimation
    y = My(:,i);  
    [~,~,res] =Spest(I,kC, kA, m, m, X, t, y,delta);
    sigCoeff=pinv(Bt'*Bt+delta*eye(size(Bt,2)))*Bt'*res.^2;
    sig2(:,i)=Bt*sigCoeff;
    Mse(i)=sqrt(mean((sig2(:,i)-truesig2).^2));
    
 end
 
 %evaluation of estimation
 ind=[mean(Mse) median(Mse) std(Mse)];
 Msig2=median(sig2,2);
 %figure comaprison of true and  estimated variance
 %plot(t,truesig2,'k-',t,mean(sig2,2),'r',t,median(sig2,2),'b')
 
%box plot comparison
% boxplot([truesig2 median(sig2,2) mean(sig2,2)])

 
end

