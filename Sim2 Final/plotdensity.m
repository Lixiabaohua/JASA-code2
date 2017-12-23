function  plotdensity(Dat, B, m, optK, optknot, delta )
% density simulation based on Monte Carlo replications and bootstrap
% sampling
xMat = Dat(:,6:7);
t = Dat(:,8);
y = Dat(:,9);
T=size(xMat,1);  
[~, Cres, Ufit] =mytest(xMat, t, y, m, optK, optknot, delta);
testM = zeros(B,1);   
testB1 = zeros(B,1);      testB2 = zeros(B,1);     testB3 = zeros(B,1);    
testB4 = zeros(B,1);      testB5 = zeros(B,1);      testB6 = zeros(B,1);  
for i = 1 : B
    % generate monte carlo replications to plot the density of test statistics
    MontY = Dat(:,1) + Dat(:, 2) .* Dat(:, 4) +  Dat(:, 3) .* Dat(:, 5) +...
                   (0.7 + (2-t)./(2+t)) .* normrnd(0,0.4,T,1);  %plot density of test statistics
    testM(i) = mytest(xMat, t, MontY, m, optK, optknot, delta);  
    z1 = Ufit + Cres .* normrnd(0,1,T,1);
    z2 = Ufit + Cres .* normrnd(0,1,T,1);
    z3 = Ufit + Cres .* normrnd(0,1,T,1);
    z4 = Ufit + Cres .* normrnd(0,1,T,1);
    z5 = Ufit + Cres .* normrnd(0,1,T,1);
    z6 = Ufit + Cres .* normrnd(0,1,T,1);
    testB1(i) = mytest(xMat, t, z1, m, optK, optknot, delta);  
    testB2(i) = mytest(xMat, t, z2, m, optK, optknot, delta);
    testB3(i) = mytest(xMat, t, z3, m, optK, optknot, delta);
    testB4(i) = mytest(xMat, t, z4, m, optK, optknot, delta);
    testB5(i) = mytest(xMat, t, z5, m, optK, optknot, delta);
    testB6(i) = mytest(xMat, t, z6, m, optK, optknot, delta);
end
% compare density of test statistics
mytest1 = [testM testB1 testB2 testB3 testB4 testB5 testB6];
xp = linspace(min(min(mytest1)),max(max(mytest1)),100);
f1 = ksdensity(mytest1(:,1),xp);
f2 = ksdensity(mytest1(:,2),xp);
f3 = ksdensity(mytest1(:,3),xp);
f4 = ksdensity(mytest1(:,4),xp);
f5 = ksdensity(mytest1(:,5),xp);
f6 = ksdensity(mytest1(:,6),xp);
f7 = ksdensity(mytest1(:,7),xp);
plot(xp, f1,'k-', xp, f2, 'r-.', xp, f3,'b-.', xp, f4, 'g-.', xp, f5, 'y-.', xp, f6, 'c-.', xp, f7, 'm-.')

end

