function  conintv = cofi_add( sampleX,myfit,myres,B,alpha,m,k,delta)
%plot  fitted curves under additive model
% myfit: fitted value of three-step spline estiamtion
% myres: residuals of three-step spline estiamtion
% sampleX: covariates,
% B: bootstrap times
% alpha: confidence level

T=size(sampleX,1);
Bc = zeros(1,B);
Bfun1 =zeros(T,B);  Bfun2 =zeros(T,B);  Bfun3 =zeros(T,B);
for i = 1:B
   %generate bootres and boot response
   randn('seed',i);
   bootres = myres .* randn(T,1);
   booty = myfit +bootres;  %bootstrap response
   [Bfun,Bconst,~,~,~]= add_est( sampleX,booty,m,k,delta );
   Bfun1(:,i) = Bfun(:,1);
   Bfun2(:,i) = Bfun(:,2);
   Bfun3(:,i) = Bfun(:,3);
   Bc(i) = Bconst;
end
[sc,~] = sort(Bc);
%i1 = ceil(B*0.5); 
i2 =ceil(B*alpha/2); i3 = ceil(B*(1-alpha/2));
conintv = [sc(i2) sc(i3)];
end