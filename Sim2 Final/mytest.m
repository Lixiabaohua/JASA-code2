function  [testValue, Cres, Ufit] = mytest(sampleX, u, y, m, Bknot,Tknot, delta)
% construct the testing statistics and centralizd error under unrestricted
% model
% sampleX: matrix of covariates
% u: rescaled time
% y: response variable
% Bknot: optimal knot for bivariate function estimation
% Tknot: optimal knot for three-step spline estimation
% delta: ridge parameters

% bivariate functions estimation
[Bfit, Bres] = BivEst(sampleX, u, y, m, Bknot(1), Bknot(2), delta);
Cres=Bres-mean(Bres); % centralized residuals
%% obtain three-step estimation 
optkC=Tknot(1); optkA = Tknot(2); optm1= Tknot(3); optN=Tknot(4);
% three-step spline estimators and error
Ufit=Spest(optN, optkC, optkA, optm1, m, sampleX, u, y, delta);
% test statistic
testValue = mean((Bfit - Ufit).^2);
end

