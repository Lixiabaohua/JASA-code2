%main procedure for NO2 air pollution data
close all ;
clear
clc;
tic
%% data and parameters
NO2Dat=xlsread('NO2.xls.xlsx', 'sheet1');
%resort observation according to time order (day number and hours of day)
sortDat = mysort(NO2Dat);
csvwrite('SNo2.csv',sortDat)
%response
y = sortDat(:, 1); 
mycorY = corr(y, sortDat(:, 2 : end));
mycorX = corr(sortDat(:, 2 : end));
SST = sum((y - mean(y)) .^ 2);
%matrix for response and covariates
xMat = sortDat(:, [2 4 5]);
%plot time series plot and autocorrelation plot for newDat
[T,p]=size(xMat); %size and dimensional of covariates
UxMat = xMat - repmat(mean(xMat), T, 1);
%cancel out the same value of covariates
rand('seed',5)
RxMat =UxMat +10^(-6) *rand(T, p);
%generate rescaled time
t=linspace(0, 1, T + 1);
t(1)=[];
t = t' ;
% optimal smoothing paramater
delta = 10 .^ -3; 
kseq = ceil(0.5 * T ^ (1/5)):ceil(2 * T ^ (1/5));
Subseq = [25 50 100];
m0seq = [2 3];  m=3;
B = 1000;
alpha = 0.05;
opt =myknot_vca( kseq, m, m0seq, Subseq, RxMat, t, y, delta ) ;
kC = opt(1); kA = opt(2); m1=opt(3);  I1 = opt(4);
%%  separablity assumption
%three-step spline estimation of varying-coefficient additive model
%three-step spline estimation of varying-coefficient additive model
Inib = StepIest(kA, m1, I1, RxMat, y, delta ) ;
[Ualp, Ubeta, Ufit, Ures, Usig] =Spest( Inib, kC, kA, m, RxMat, t, y, delta);
Ur = 1 - sum(Ures .^ 2)/SST;
% bivariate functions estimation
optBK = Bivknot(kseq, m, RxMat, t, y, delta);  
[Bfit, Bres] = BivEst(RxMat, t, y, m, optBK(1), optBK(2), delta);
[testV, Cres] = mytest(Ufit, RxMat, t, y, m, optBK, delta);
[pv, cri] = bootCri(B, testV, Ufit, Cres, t, RxMat, m, optBK, alpha, delta); % p value
%% model identification
bound = 10^(-2);    
lamseq =0 : 0.05 : 1;
museq =0 : 0.05:1 ;
[lam, mu, varyInd, linearInd, Pnalp, Pnbeta] =...
     optimtune(Ubeta, RxMat, t, y, lamseq, museq,  bound, delta, p, kC, kA, m);
%conclusion: reduced model 
[~, ~, Rfit, Rres, Rsig] = RM( RxMat, t, y, m, kC, kA, delta );
Rr = 1-sum(Rres .^ 2)/SST;
%plot under VCAM
B = 1000;  %bootstrap times
alpha =0.05;   %confidence level 
%constant estimation and confidence bands under additive model
[~, const, Afit, Ares] = add_est(RxMat, y, m, kA, delta );
conintv = cofi_add(RxMat, Afit, Ares, B, alpha, m, kA, delta) ;
%plot for reduced model: Figure 4 in Example 3
plot_RM(RxMat, t, Rfit, Rres, B, alpha, m, kC, kA, delta, const, conintv)
n = 30;
[err, yA] = RM_test( T, RxMat, t, y, n, m, kC, kA, delta);
preMse = sqrt(mean(err .^ 2));
pt = 471:500;
% plot one-step forword forecasting
plotpre(pt, Rfit, yA, err)
mytime = tic;