function plotRMest( sampleX, u, myfit, myres, B, alpha, m1, m, kC, kA, I1, n1, I2, n2,  delta)
          
%plot  fitted curves under additive model
% myfit: fitted value of three-step spline estiamtion
% myres: residuals of three-step spline estiamtion
% sampleX: covariates,
% B: bootstrap times
% alpha: confidence level
tseq = linspace( quantile(u, 0.01), quantile( u, 0.99), 30);
x1seq = linspace( quantile( sampleX(:,1), 0.01), quantile( sampleX(:,1), 0.99 ), 30) ;
x2seq = linspace( quantile( sampleX(:,2), 0.01),   quantile( sampleX(:,2), 0.99), 30) ;
T=size(sampleX, 1);  Lt = length(tseq); Lx=length(x1seq);
Balp0 = zeros(Lt, B);   Balp1 = zeros(Lt, B);  
Bbeta1 =zeros(Lx, B);  Bbeta2 =zeros(Lx, B);  
Cres = myres - mean(myres) ;
for i = 1:B
   booty = myfit + Cres .* normrnd(0, 1, T, 1) ;  %bootstrap response
   Inib = StepIest(kA, m1, I1, n1, I2, n2, sampleX, booty, delta ) ;
   [ ~, ~, coefC, coefA, a, c] = RM( Inib, kC, kA, m, sampleX, u, booty, delta ) ;
   [Bhalp, Bbeta] = pred_RMest(x1seq', x2seq', tseq', sampleX, u, coefC, coefA, a, c, m, kC, kA ) ;
   Balp0(:,i)  =    Bhalp(:, 1) ;
   Balp1(:,i)  =    Bhalp(:, 2);
   Bbeta1(:,i)  =  Bbeta(:, 1);
   Bbeta2(:,i)  =  Bbeta(:, 2);
end
[salp0,~] = sort(Balp0,2);
[salp1,~] = sort(Balp1,2) ;
[sbeta1,~] = sort(Bbeta1,2);
[sbeta2,~] = sort(Bbeta2,2) ;
i1 = ceil(B*0.5); i2 =ceil(B*alpha/2); i3 = ceil(B*(1-alpha/2));
efun1 = [ mean(sbeta1,2) sbeta1(:,i1) sbeta1(:,i2) sbeta1(:,i3)];
efun2 = [ mean(sbeta2,2) sbeta2(:,i1) sbeta2(:,i2) sbeta2(:,i3)];


plot(tseq, salp0(:,i1), 'r-.',  tseq, salp0(:,i2),'b--', tseq, salp0(:,i3), 'b--', 'LineWidth', 1)
xlabel('Rescaled Time')
ylabel('\alpha_{0}')
xlim([-0.05 1.05])


plot(tseq, salp1(:, i1), 'r-.',  tseq, salp1(:, i2), 'b--', tseq, salp1(:, i3), 'b--', 'LineWidth', 1)
xlabel('Rescaled Time')
ylabel('\alpha_{1}')
xlim([-0.05 1.05])




plot(x1seq, efun1(:, 2),'r-.',...
       x1seq, efun1(:, 3),'b--', x1seq, efun1(:, 4),'b--', 'LineWidth',1)
xlabel('Y_{t-1}')
ylabel('\beta_{1}')
xlim([-0.02 0.06])


plot(x2seq, efun2(:, 2),'r-.',...
       x2seq,efun2(:, 3), 'b--', x2seq, efun2(:, 4), 'b--', 'LineWidth', 1)
xlabel('R_{t-1}')
ylabel('\beta_{2}')
xlim([-0.06 0.06])

end




