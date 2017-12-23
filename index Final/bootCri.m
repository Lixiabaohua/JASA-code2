function  [pv, cri] = bootCri(B1, testV, Ufit, Cres, t, xMat, m, optK, delta, alpha)
% generate wild bootstrap sampling
%Q: size of bootstrap resampling
%Dat: DGP
% m and m1seq are order of B-spline in Step II (III) estimation and Step I
% estimation
% Nseq: sequence of segment length in Step I estimation
% delta: ridge parameter
% B1: size of bootstrap sampling to decide critical value
% optK for bivariate estimation; optknot for VCAM
T= length(t) ;
%%  compute p value
testBV = zeros(1,B1); 
p1=zeros(1,B1);   
for i = 1 : B1
    booty     = Ufit  + Cres .* normrnd(0, 1, T, 1);    % bootstrap response
    testBV(i) = mytest(Ufit,  xMat, t,  booty, m, optK, delta) ; % value of test statistic
    p1(i)       = p1(i) + (testBV(i) > testV);             % frequency of rejecting H0    
end
stest = sort(testBV);
cri = stest(B1*(1-alpha));
pv = mean(p1); 
end
  
