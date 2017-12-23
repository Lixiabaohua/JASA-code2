  function [ outMat, Sy] = mydata3(T,Q)
%generate the data with sample size T
%varying-coefficient function
%% varying-coefficient function
t=linspace(0,1,T+1);
t(1)=[];
t = t' ;
alpha0=1.5*t+2*cos(2*pi*t);
alpha1 = 1.2*(1-t).*cos(2*pi*t) +1;
%alpha1 = 1.3*t.*sin(2*pi*t) +1;
alpha2 =  1.5*(1-t).^2;
alpha3 = ones(T,1) ;
alpha4 = 0.5*sin(2*pi*t)+1;
%standardized
alp = [alpha0 alpha1 alpha2 alpha3 alpha4];
aMat = alp./repmat([1 mean(alp(:,2:end))],T,1) ;
%% generatr independently local stationary covariates
x1 = autoreg1( T,0.6,0.4 ) ;
x2 = autoreg1( T,0.4,0.4) ;
x3 = autoreg1( T, 0.5, 0.4) ;
x4 = autoreg1(T,  0.7,0.6);
x1 = x1- mean(x1);
x2 = x2- mean(x2);
x3 = x3 -mean(x3);
x4 = x4- mean(x4);
%% additive function
beta1 =  0.7*sin(pi * x1/2)-0.5 * x1 .* (2-x1).^2;
beta2 = 3 * cos(pi * x2/2).*x2;
beta3 = 2 * x3 .* (1+x3) ;
beta4 = x4;
%centering
beta1 = beta1 - mean(beta1);
beta2 = beta2 - mean(beta2);
beta3 = beta3 - mean(beta3);
bMat = [beta1 beta2 beta3 beta4];
%% response variable
y=aMat(:,1)+aMat(:,2).*bMat(:,1)+aMat(:,3).*bMat(:,2)+...
    aMat(:,4).*bMat(:,3) + aMat(:,5).*bMat(:,4)+normrnd(0,0.5,T,1);
outMat = [aMat bMat x1 x2 x3 x4 t y]; %15 columns
%% generate Monte Carlo replications
Sy = zeros(T,Q);
for i = 1: Q
    Sy(:,i) = aMat(:,1)+aMat(:,2).*bMat(:,1)+aMat(:,3).*bMat(:,2)+...
                     aMat(:,4).*bMat(:,3) + aMat(:,5).*bMat(:,4)+normrnd(0,0.5,T,1);
end
end

