%% Code for Example 2: multiplicity testing
close all
clear 
clc
tic
%% parameters
% size of bootstrap samples, former for criterion value and p-value; the
% latter for power
B1 = 1000;  B2 =600; 
delta = 0.001;
% level
alpha = 0.05;  
thet = [0  0.5  1.0 1.5 2 2.5];
m = 3; m1seq = [2 3]; Nseq = [25 50 100 ];
%% load data
% Dat1
load('Dat1');  load('normR11'); load('normR12'); 
power1 = myPower(B1, B2, Dat1, normR11, normR12, m, m1seq, Nseq, delta, thet, alpha);
% Dat 2
load('Dat2');  load('normR21'); load('normR22'); 
power2 = myPower(B1, B2, Dat2, normR21, normR22, m, m1seq, Nseq, delta, thet, alpha);
% Dat 4
load('Dat4');  load('normR41'); load('normR42'); 
power3 = myPower(B1, B2, Dat4, normR41, normR42, m, m1seq, Nseq, delta, thet, alpha);
%% Print 
% print power Table 
MatPower = [power1' power2' power3'];
% plot  power function
plot(thet,power1,'k:',thet,power2,'b--',thet,power3,'r-.','LineWidth',1)
xlim([-0.1 2.6])
ylim([-0.1 1.1])
legend('T=200','T=300','T=500','Location', 'southeast')
% plot density
load('Dat3')
[~, ~, optK, optknot] = myPower(B1, B2, Dat3, normR31, normR32, m, m1seq, Nseq, delta, thet, alpha);
plotdensity(Dat3, 1000, m, optK, optknot, delta)
mytime = toc