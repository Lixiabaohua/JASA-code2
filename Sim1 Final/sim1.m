 %% code for Example 1 in paper
 % Table 1 in the body: comparisons of estimations
 % Table 2 (Table 1  in the supplementary materials) : comparisons of predictions
 % Table 3 (Table 2  in the supplementary materials) : comparisons of multi-step iterations
 % Figure 1 (in the body) : pointwise confidence bands and simultaneous confidence bands
 % Figure 2 (in the body) : bootstrap based confidence bands
 % Figure 3 (Figure 1 in the supplementary materials): histgram of  bootstrap statistics
 % Figure 4 (Figure 2 in the supplementary materials): comaprisons of estimated surfaces 
close all
clear
clc
tic
%% parameters
%optional cases in each group
Subseq = [25 50 100 150]; 
% order of B-spline in Step II and Step III estimation
m = 3;  m1seq = [2 3];
delta = 10^(-3);
Q =1000;
%% Dat 1
load('Dat1')
load('My1')
[ MseU1, MseO1, MseA1, MseV1, Msey1, knot1 ] = Msecomp(Q, Dat1, My1, m, m1seq, Subseq, delta) ; 
[MseI1,MseyI1] = IterMse( Q, Dat1, My1, m, m1seq, Subseq, delta);
%results for Table 1
Res1 = [MseU1' MseO1' [MseV1(:,1:3) [nan nan; nan nan] MseV1(:,4)]'  ...
          [MseA1(:,1) [nan nan; nan nan] MseA1(:,2:4)]' ];
%results for Table 2    
outMat1= mse_Pre(Q, Dat1, My1, 5,  knot1(1), knot1(2), knot1(3), m, knot1(4), delta ) ;
outMat2= mse_Pre(Q, Dat1, My1, 10,  knot1(1), knot1(2), knot1(3), m, knot1(4), delta ) ;
%results for Table 3 
Iter1=[MseI1; MseyI1];
%% Dat 2
load('Dat2')
load('My2')
[ MseU2, MseO2, MseA2, MseV2, Msey2, knot2 ] = Msecomp(Q, Dat2, My2, m, m1seq, Subseq, delta) ; 
[MseI2, MseyI2] = IterMse( Q, Dat2, My2, m, m1seq, Subseq, delta);
% results for Table 1
Res2 = [MseU2' MseO2' [MseV2(:,1:3) [nan nan; nan nan] MseV2(:,4)]'  ...
            [MseA2(:,1) [nan nan; nan nan] MseA2(:,2:4)]' ];
% results for Table 2        
outMat3= mse_Pre(Q, Dat2, My2, 5,  knot2(1), knot2(2), knot2(3), m, knot2(4), delta ) ;
outMat4= mse_Pre(Q, Dat2, My2, 10,  knot2(1), knot2(2), knot2(3), m, knot2(4), delta ) ;
%results for Table 3 
Iter2=[MseI2; MseyI2];
%% Dat 3
load('Dat3')
load('My3')
[ MseU3, MseO3, MseA3, MseV3, Msey3, knot3 ] = Msecomp(Q, Dat3, My3, m, m1seq, Subseq, delta) ; 
[MseI3, MseyI3] = IterMse( Q, Dat3, My3, m, m1seq, Subseq, delta);
% results for Table 1
Res3 = [MseU3' MseO3' [MseV3(:,1:3) [nan nan; nan nan] MseV3(:,4)]'  ...
          [MseA3(:,1) [nan nan; nan nan] MseA3(:,2:4)]' ];
% results for Table 2
outMat5= mse_Pre(Q,Dat3, My3, 5,  knot3(1), knot3(2), knot3(3), m, knot3(4), delta ) ;
outMat6= mse_Pre(Q, Dat3, My3, 10,  knot3(1), knot3(2), knot3(3), m, knot3(4), delta ) ;
%results for Table 3 
Iter3=[MseI3; MseyI3];
%% Final printout 
%Table 1
MseMat = [Res1; Res2; Res3];
%print optimal parameter in Table 1 
paraMat = [knot1; knot2; knot3];
%Table 2 
PreMat = [outMat1;  outMat2; outMat3;  outMat4; outMat5;   outMat6];                
%Table 3 
IterMat=[Iter1; Iter2; Iter3];
%% Plot 
load('Dat4')
% sample size
T=size(Dat4, 1);
kseq = ceil(0.5 * T ^ (1/5)) : ceil(2 * T ^ (1/5));
myplot( Dat4, para4, T, Subseq, kseq, m, m1seq, Q, delta )
mytime = toc



