function  [V_power, pv, optK, optknot] = myPower(B1, B2, Dat, normR1, normR2, m, m1seq, Nseq, delta, thet, alp)
% generate wild bootstrap sampling
%Q: size of bootstrap resampling
%Dat: DGP
% m and m1seq are order of B-spline in Step II (III) estimation and Step I
% estimation
% Nseq: sequence of segment length in Step I estimation
% delta: ridge parameter
% B1: size of bootstrap sampling to decide critical value
% B2: size of bootstrap sampling to decide power(rejection frequency)
%% dat preparation
xMat = Dat(:,6:7);
t = Dat(:,8);
y = Dat(:,9);
T=size(xMat,1);  
%% decide optimal parameters
kseq = ceil(0.5*T^(1/5)) : ceil(2*T^(1/5));
% obtain optimal knots of bivariate spline estimation 
optK = Bivknot(kseq, m, xMat, t, y, delta);
% obtain optimal knots of three-step spline estimation 
optknot = myknot_vca(kseq, m, m1seq, Nseq, xMat, t, y, delta) ;
%%  decide the critical based on B1 bootstrap sampling
[testV, Cres, Ufit] = mytest(xMat, t, y, m, optK, optknot, delta);
testBV = zeros(1,B1); 
p1=zeros(1,B1);
testM      =  zeros(1,B2);  % restore the values of test statistics based on Monte Carlo replications
for i = 1:B1
    booty     = Ufit  + Cres .* normR1(:,i);     % bootstrap response
    testBV(i) = mytest(xMat, t, booty, m, optK, optknot, delta); % value of test statistic
    p1(i)       = p1(i) + (testBV(i) > testV);             % frequency of rejecting H0
     % generate monte carlo replications to plot the density of test statistics
    MontY = Dat(:,1) + Dat(:, 2) .* Dat(:, 4) +  Dat(:, 3) .* Dat(:, 5) +...
                   (0.7 + (2-t)./(2+t)) .* normrnd(0,0.4,T,1);  %plot density of test statistics
    testM(i) = mytest(xMat, t, MontY, m, optK, optknot, delta);  % value of test statistic
    
end
stest = sort(testBV);
cri = stest(B1*(1-alp)); % criterion value based on 1000 bootstrap samples
pv = mean(p1);
%% compute power for nine alternative hypothesis

% restore the frequency of rejection
s1 = 0;  s2 = 0;   s3 = 0;  s4 =0;   s5 =0;   s6 = 0;    %s7 = 0;   s8 = 0;  s9 = 0; s10=0;
% restore the values of test statistics based on bootstrap sampling
for i = 1 : B2
    % bootstrap response under different hypothesting
    y1 =  Ufit + Cres .* normR2(:,i) + thet(1) * (1.2*xMat(:, 1) + 0.8*t);
    y2 =  Ufit + Cres .* normR2(:,i) + thet(2) * (1.2*xMat(:, 1) + 0.8*t);
    y3 =  Ufit + Cres .* normR2(:,i) + thet(3) * (1.2*xMat(:, 1) + 0.8*t);
    y4 =  Ufit + Cres .* normR2(:,i) + thet(4) * (1.2*xMat(:, 1) + 0.8*t);
    y5 =  Ufit + Cres .* normR2(:,i) + thet(5) * (1.2*xMat(:, 1) + 0.8*t);
    y6 =  Ufit + Cres .* normR2(:,i) + thet(6) * (1.2*xMat(:, 1) + 0.8*t);
    testBV1 = mytest(xMat, t, y1, m, optK, optknot, delta); 
    testBV2 = mytest(xMat, t, y2, m, optK, optknot, delta); 
    testBV3  = mytest(xMat, t, y3, m, optK, optknot, delta); 
    testBV4     = mytest(xMat, t, y4, m, optK, optknot, delta); 
    testBV5  = mytest(xMat, t, y5, m, optK, optknot, delta); 
    testBV6     = mytest(xMat, t, y6, m, optK, optknot, delta); 
    s1 = s1 + (testBV1   > cri);
    s2 = s2+  (testBV2  > cri);
    s3 = s3 + (testBV3 > cri);
    s4 = s4+  (testBV4  > cri);
    s5 = s5+  (testBV5  > cri);
    s6 = s6 + (testBV6  > cri);
end
V_power =[s1 s2 s3 s4 s5 s6 ]/B2;

end
  
