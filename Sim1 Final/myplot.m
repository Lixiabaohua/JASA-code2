function  myplot( Dat5,para5,T,Subseq,kseq,m,m1seq,Q,delta )
% plot Figures of simulation 1
% choose optimal knots for varying-coefficient additive model
[optknot,~] = myknot_vca(kseq,m,m1seq,Subseq,Dat5(:,6:7),Dat5(:,8),Dat5(:,9),delta) ;
kC=optknot(1); kA = optknot(2); optm1= optknot(3); optN=optknot(4);
% three-step spline estimators
[ealp,ebet,res] =Spest( optN,kC, kA, optm1, m,Dat5(:,6:7),Dat5(:,8),Dat5(:,9),delta);
Cres = res - mean(res);
% generate variance of random errors
Bt= rspline(Dat5(:,8), Dat5(:,8) ,m, kC); 
sigCoeff=pinv(Bt'*Bt+delta*eye(size(Bt,2)))*Bt'*res.^2;
sig2=Bt*sigCoeff;
%plot Figure 1: bands
%kA=3;
plotbandC(Dat5, ealp, ebet, sig2, Cres, kC, kA, m, optm1, optN, Q, delta )
% plot Figure 2: histogram and QQ plot
plotboot( Dat5, ealp, ebet, sig2, T, Q, optN, kC, kA, optm1, m, delta)
% plot Figure 3: surface plot
plotsurf(Dat5,para5)





end

