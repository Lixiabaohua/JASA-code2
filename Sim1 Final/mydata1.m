function [outMat,stdpara] = mydata1(T)
%generate the data with sample size T
%rescaled time
t=linspace(0,1,T+1);
t(1)=[];
t = t' ;
%% varying-coefficient functions
alpha0 = 1.5*cos(pi*t*2);
alpha1 = 1.3*t.*sin(2*pi*t) +1;
alpha2 =2*sin(1.5*pi*t)-1.2*(t-0.5).*(1-t)+1;
%regularization
norm1 =  mean(alpha1);
norm2 =  mean(alpha2);
alp1 = alpha1/norm1;
alp2 = alpha2/norm2;
aMat = [alpha0 alp1 alp2];
%% additive functions
%generatr independently local stationary covariates
x1 = autoreg1( T,0.6,0.5 ) ;
x2 = autoreg1( T,0.8,0.3) ;
%additive function
beta1 = 3*sin(pi*x1/2) - (1-x1).*x1;
beta2 = 2*cos(pi*x2/2) +1.8*sin(pi*x2/3);
%centering
c1 = mean(beta1);
c2 = mean(beta2);
beta1 = beta1 - c1;
beta2 = beta2 - c2;
%% response variable
y=aMat(:,1)+aMat(:,2).*beta1+aMat(:,3).*beta2+(0.7+(2-t)./(2+t)).*normrnd(0,0.4,T,1);
outMat = [aMat beta1 beta2  x1 x2 t y];
stdpara = [norm1 norm2 c1 c2];
end

