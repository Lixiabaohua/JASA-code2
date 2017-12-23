% data between 2000/1/1 anf 2011/12/30
%
close all; 
clear; 
clc;
tic
%% data and parameters
%read oroiginal data
Data = csvread('ixica.csv', 1, 8);
% plot time series (Figure 4)
MyDat = Data(1 : 3019,:);
[T,p]=size(MyDat);
plotseries(MyDat ) 
y = Data(2 : T,2);        %response
x1 = Data(1 : T - 1,2);
x2 = Data(1 : T - 1,1)/100;  %covariate unit :percentile%
x = [x1 x2];
T=length(y);
SST = sum((y - mean(y)).^2);
%generate rescaled time
t=linspace(0, 1, T + 1);
t(1)=[];
t = t' ;
%centering covariates on [0,1]
ux = x - repmat(mean(x), T, 1);
%cancel out the same value of covariates
rand('seed', 5)
RxMat =ux + 10^(-6) * rand(T, p);
kseq = ceil(0.5 * T ^ (1/5)):ceil(2 * T ^ (1/5)); 
m0seq = [2 3]; 
delta = 10^(-2);
m=3;
B1=1000;  
alpha=0.05;
%choose optimal smoothing parameter
opt =myknot_vca1( kseq, m, m0seq, RxMat, t, y, delta ) ;
kC=opt(1); kA =opt(2); m1=opt(3);  I1 = opt(4); n1 =opt(5); I2=opt(6); n2=opt(7);
%% separablity assumption
%three-step spline estimation of varying-coefficient additive model
Inib = StepIest(kC, m1, I1, n1, I2, n2, RxMat, y,delta ) ;
[Ualp,Ubeta, Ufit, Ures, Usig] =Spest( Inib, kC, kA, m, RxMat, t, y, delta);
Ur = 1-sum(Ures .^ 2)/SST; 
optBK = Bivknot(kseq, m, RxMat, t, y, delta);  
[testV, Cres] = mytest(Ufit, RxMat, t, y, m, optBK, delta);
[pv, cri] = bootCri(B1, testV, Ufit, Cres, t, RxMat, m, optBK, delta, alpha);
[Bfit, Bres] = BivEst(RxMat, t, y, m, optBK(1), optBK(2), delta);
%% model comparison
%Varying-coefficient model
%optimal knots for varying-coefficient model
k_vc = knot_vc( kseq, m, RxMat, t, y, delta ) ;
[ ~,res_vc,sig ] =vcm(RxMat, t, y, m, k_vc, delta);
R_vc = 1-sum(res_vc.^2)/SST;
% Additive model
k_add = knot_add(kseq, RxMat, y, m, delta );
[~,~,~,res_add]=add_est( RxMat, y, m, k_add, delta );
R_add = 1-sum(res_add.^2)/SST;
%% model identification
bound = 10^(-3);    
lamseq =0 : 0.01 : 1;
museq  = 0 : 0.01 : 1 ;
[lam, mu, varyInd, linearInd, Pnalp, Pnbeta] =...
     optimtune(Ubeta, RxMat, t, y, lamseq, museq, bound, delta, p, kC, kA, m);
% for reduced model
[Rmres,Rmfit] = RM( Inib, kC, kA, m, RxMat, t, y, delta ) ;

df = 2 * (m + kC) + 2 * (m + kA);
Rmsig = sum(Rmres .^ 2)/(T-df);
Rr2 = 1-sum(Rmres .^ 2)/SST;
%% plot
% plot for reduced model
plotRMsurf(Inib, RxMat, y, t, kC, kA, m, delta)
plotRMest( RxMat, t, Ufit, Ures, B1, alpha, m1, m, kC, kA, I1, n1, I2, n2,  delta)          
%% prediction
n = 50;   
[err, ypre] = Sp_test(T, RxMat, t, y, n, kC, kA, m1, m, delta);
preMse1 = sqrt(mean(err.^2));
% plot one-step forescasting (Figure 7)
plotpre(Rmfit, ypre, err,  T)
mytime = toc