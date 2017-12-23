function plotboot( Dat,ealp,ebet,eSig,T,Q,optN,kC,kA,optm1,m,delta)
% construct wide bootstrap sampling
% ealp and ebet are three-step spline estimators
% eSig: estimated standard deviation of random error
%T,kC,kA are sample size, and optimal knots number
%Dat are sample(9 column)
%Q: size of bootstrap sampling
% m and optm1 are order of B-spline
%% data
x=Dat(:,6:7);
u=Dat(:,8);
%% bootstrap confidence bands
Balp0=zeros(T,Q);   Balp1=zeros(T,Q); Balp2=zeros(T,Q);
Bbet1=zeros(T,Q);   Bbet2=zeros(T,Q);
for i=1:Q
   y= ealp(:,1)+ealp(:,2).*ebet(:,1)+ealp(:,3).*ebet(:,2)+eSig.*normrnd(0,1,T,1);
   [halp,hbeta] =   Spest( optN,kC, kA, optm1, m,x,u,y,delta);
   Balp0(:,i)=halp(:,1);
   Balp1(:,i)=halp(:,2);
   Balp2(:,i)=halp(:,3);
   Bbet1(:,i)=hbeta(:,1);
   Bbet2(:,i)=hbeta(:,2);
end
%%  different distribution plot
[~,ind]=sort(x);
% plot  histgram at three 
%histfit(Balp0(T*0.25,:),[],'kernel')
%title('Histogram and Kernel Density of Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

histfit(Balp0(T*0.5,:), [],'kernel')
title('Histogram and Kernel Density of Bootstrap Statistic')
xlabel('Bootstrap Statistic')
ylabel('Frequency')

%histfit(Balp0(T*0.75,:), [],'kernel')
%title('Histogram and Kernel Density of Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

%histfit(Balp1(T*0.25,:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

histfit(Balp1(T*0.5, :), [],'kernel')
title('Histogram and Kernel Density of Bootstrap Statistic')
xlabel('Bootstrap Statistic')
ylabel('Frequency')

%histfit(Balp1(T*0.75,:), [],'kernel')
%title('Histogram and Kernel Density of Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

%histfit(Balp2(T*0.25,:), [],'kernel')
%title('Histogram and Kernel Density of Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

histfit(Balp2(T*0.5,:), [],'kernel')
title('Histogram and Kernel Density of  Bootstrap Statistic')
xlabel('Bootstrap Statistic')
ylabel('Frequency')

%histfit(Balp2(T*0.75,:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

%histfit(Bbet1(ind(T*0.25,1),:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

histfit(Bbet1(ind(T*0.5,1),:), [],'kernel')
title('Histogram and Kernel Density of  Bootstrap Statistic')
xlabel('Bootstrap Statistic')
ylabel('Frequency')

%histfit(Bbet1(ind(T*0.75,1),:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

%histfit(Bbet2(ind(T*0.25,2),:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

histfit(Bbet2(ind(T*0.5,2),:), [],'kernel')
title('Histogram and Kernel Density of  Bootstrap Statistic')
xlabel('Bootstrap Statistic')
ylabel('Frequency')

%histfit(Bbet2(ind(T*0.75,2),:), [],'kernel')
%title('Histogram and Kernel Density of  Bootstrap Statistic')
%xlabel('Bootstrap Statistic')
%ylabel('Frequency')

% probability plot 
%normplot(Balp0(T*0.25,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp0(T*0.5,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp0(T*0.75,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp1(T*0.25,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp1(T*0.5,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp1(T*0.75,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp2(T*0.25,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp2(T*0.5,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Balp2(T*0.75,:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet1(ind(T*0.25,1),:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet1(ind(T*0.5,1),:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet1(ind(T*0.75,1),:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet2(ind(T*0.25,2),:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet2(ind(T*0.5,2),:))
%xlabel('Bootstrap Quantiles')
%box on

%normplot(Bbet2(ind(T*0.75,2),:))
%xlabel('Bootstrap Quantiles')
%box on

%  Q-Q plot
%qqplot(Balp0(T*0.25,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

qqplot(Balp0(T*0.5,:))
title('QQ Plot of Bootstrap Statistic versus Standard Normal')
ylabel('Bootstrap Quantiles')
box on

%qqplot(Balp0(T*0.75,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

%qqplot(Balp1(T*0.25,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

qqplot(Balp1(T*0.5,:))
ylabel('Bootstrap Quantiles')
title('QQ Plot of Bootstrap Statistic versus Standard Normal')
box on

%qqplot(Balp1(T*0.75,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

%qqplot(Balp2(T*0.25,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

qqplot(Balp2(T*0.5,:))
ylabel('Bootstrap Quantiles')
title('QQ Plot of Bootstrap Statistic versus Standard Normal')
box on

%qqplot(Balp2(T*0.75,:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

%qqplot(Bbet1(ind(T*0.25,1),:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

qqplot(Bbet1(ind(T*0.5,1),:))
ylabel('Bootstrap Quantiles')
title('QQ Plot of Bootstrap Statistic versus Standard Normal')
box on

%qqplot(Bbet1(ind(T*0.75,1),:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

%qqplot(Bbet2(ind(T*0.25,2),:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on

qqplot(Bbet2(ind(T*0.5,2),:))
ylabel('Bootstrap Quantiles')
title('QQ Plot of Bootstrap Statistic versus Standard Normal')
box on

%qqplot(Bbet2(ind(T*0.75,2),:))
%ylabel('Bootstrap Quantiles')
%title('QQ Plot of Bootstrap Statistic versus Standard Normal')
%box on
end

