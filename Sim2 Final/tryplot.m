% plot density
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
load('Dat3')
[~, ~, optK, optknot] = myPower(B1, B2, Dat3, normR31, normR32, m, m1seq, Nseq, delta, thet, alpha);
plotdensity(Dat3, 1000, m, optK, optknot, delta)