function plot_RM( sampleX, u, myfit, myres, B, alpha, m, kC, kA, delta, c, cint)
%plot  fitted curves under additive model
% myfit: fitted value of three-step spline estiamtion
% myres: residuals of three-step spline estiamtion
% sampleX: covariates,
% B: bootstrap times
% alpha: confidence level
tseq = linspace( quantile(u, 0), quantile( u, 1), 30);
x1seq = linspace( quantile( sampleX(:,1), 0.025), quantile( sampleX(:,1), 0.975 ), 30) ;
x2seq = linspace( quantile( sampleX(:,2), 0.025),   quantile( sampleX(:,2), 0.975), 30) ;
x3seq = linspace( quantile( sampleX(:,3), 0.025),   quantile( sampleX(:,3), 0.975), 30) ;
T=size(sampleX, 1);  Lt = length(tseq); Lx=length(x1seq); 
Bc = zeros(Lt ,B);
Bfun1 =zeros(Lx, B);  Bfun2 =zeros(Lx, B);  Bfun3 =zeros(Lx, B);
Cres = myres - mean(myres);
for i = 1 : B
   %generate bootres and boot response
   booty = myfit +Cres .* normrnd(0, 1, T,1);  %bootstrap response
   [~, ~,~, ~, ~, coeff, mymean, up0]= RM( sampleX, u, booty, m, kC, kA, delta );
   [Balp, Bbeta] = pred_RMest(x1seq', x2seq', x3seq', tseq', sampleX, u, coeff, mymean,up0,  m, kC, kA ) ;
   Bfun1(:,i) = Bbeta(:,1);
   Bfun2(:,i) = Bbeta(:,2);   
   Bfun3(:,i) = Bbeta(:,3);
   Bc(:,i) = Balp;
end
[salp0,~] = sort(Bc,2);
[sbeta1,~] = sort(Bfun1,2);
[sbeta2,~] = sort(Bfun2,2) ;
[sbeta3,~] = sort(Bfun3,2) ;
i1 = ceil(B*0.5); i2 =ceil(B*alpha/2); i3 = ceil(B*(1-alpha/2));
efun1 = [ mean(sbeta1,2) sbeta1(:, i1) sbeta1(:, i2) sbeta1(:, i3)];
efun2 = [ mean(sbeta2,2) sbeta2(:, i1) sbeta2(:, i2) sbeta2(:, i3)];
efun3 = [ mean(sbeta3,2) sbeta3(:, i1) sbeta3(:, i2) sbeta3(:, i3)];

plot(tseq, salp0(:, i1), 'r-.',  tseq, salp0(:,i2), 'b--', tseq, salp0(:,i3), 'b--', 'LineWidth', 1)
xlabel('rescaled time')
ylabel('\alpha_{0} for constant')
xlim([-0.05 1.05])
hold on 
plot([-0.05 1.05],[c c],'k-','LineWidth',1)
plot([-0.05 1.05],[cint(1) cint(1)],'k:','LineWidth',1)
plot([-0.05 1.05],[cint(2) cint(2)],'k:','LineWidth',1)
hold off

plot(x1seq, efun1(:, 2),'r-.',...
       x1seq, efun1(:, 3),'b--', x1seq, efun1(:, 4), 'b--', 'LineWidth', 1)
xlabel('log car numbers per hour')
ylabel('\beta_{1}')



plot(x2seq, efun2(:, 2), 'r-.',...
       x2seq, efun2(:, 3), 'b--',  x2seq, efun2(:, 4), 'b--',  'LineWidth', 1)
xlabel('wind speed')
ylabel('\beta_{2}')
xlim([-3 4.2])


%subplot(3,2,4)
plot(x3seq, efun3(:, 2), 'r-.',...
       x3seq, efun3(:, 3), 'b--', x3seq, efun3(:, 4), 'b--',  'LineWidth', 1)
xlabel('difference of temperature')
ylabel('\beta_{3}')
xlim([-3 2.2])


%subplot(3,2,5)
plot(myres,'b.')
hold on
plot([0 500],[0 0],'r--','LineWidth',1)
hold off

%subplot(3,2,6)
qqplot(myres)
box on 
end

