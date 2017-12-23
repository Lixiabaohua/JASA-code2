%% Code for Example 3: two-stage model identification procedure
% print out 
% Table 3 in the body: comparisons of estimators
% Table 4 in the body: performance of model identification
close all
clear
clc
tic
%% Parameters
%optional cases in each group1
Subseq =   [25 50 100 150];
Subseq1 = [ 40 50 80 100]; 
% order of B-spline in Step II and Step III estimation
m = 3;  m1seq = [2 3]; 
delta = 10^(-3); bound = 10.^-3;
Q = 1000;
c1 = 0.6; c2 = 1.2; c0 = 1.15;
%% Dat 1 
load('Dat1')
load('My1')
[MseU1, MseOV1, MseOA1, MseP1,  tune1, test1, bias1] = ...
              Msecomp( Q, Dat1, My1, m, m1seq, Subseq, delta, bound, c1, c0);                                                                                                                                                                                                                                                                                                                                                                             delta = 10^(-3); gnum = 25;
comp1 = [MseU1; MseP1; [MseOV1(:,1 : 3) [nan; nan] MseOV1(:,4)] MseOA1 [nan; nan]];  
%% Dat 2
load('Dat2')
load('My2')
[MseU2, MseOV2, MseOA2,MseP2 , tune2, test2, bias2] = ...
             Msecomp(Q, Dat2, My2, m, m1seq, Subseq, delta, bound, c1, c2);                                                                                                                                                                                                                                                                                                                                                                     delta = 10^(-3); gnum = 25;
comp2 = [MseU2; MseP2; [MseOV2(:,1  : 3) [nan; nan] MseOV2(:,4)] MseOA2 [nan; nan]];
%% Dat 3
load('Dat3')
load('My3')
[MseU3, MseOV3, MseOA3, MseP3 , tune3, test3, bias3] = ...
     Msecomp( Q, Dat3, My3, m, m1seq, Subseq, delta, bound, c1, c2);
comp3 = [MseU3; MseP3; [MseOV3(:, 1 : 3) [nan; nan] MseOV3( :, 4) ] MseOA3 [nan; nan]]; 
%% Print
%print Table 3 in Example 3
mseMat = [comp1'; comp2'; comp3'];
% print Table 4  in Example 3
testMat = [test1; test2; test3 ];
proTest = testMat/10;
%optimal parameters under data sets
paramat = [tune1; tune2; tune3];
biasmat = bias3;
mytime = toc
