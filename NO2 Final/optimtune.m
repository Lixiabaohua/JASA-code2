function [optlam, optmu, varyInd,linearInd,Pnalp,Pnbeta] = ...
    optimtune( funbeta,sampleX,u,y,lamseq,museq,bound,delta,p,k1,k2,m)
%choose optimal tuning parameter
%k1: knots for varying-coefficient 
%k2: knots for additive 
%Stage I
mybic1 = zeros(1,length(lamseq));
for i = 1:length(lamseq)
    [~,~,mybic1(i)]= tune_StageI(funbeta,k1,m, lamseq(i),u,y,bound,p,delta);
end
[~,s1ind] =sort(mybic1);
optlam= lamseq(s1ind(1));
%norm of penalized varying-coefficient function estimation under optimal tuning parameter
[Palp,Pnalp]= tune_StageI(funbeta,k1,m, optlam,u,y,bound,p,delta);
varyInd = find(Pnalp<bound);
%Stage II 
mybic2= zeros(1,length(museq));
for i = 1:length(museq)
[~,~,~,mybic2(i)]= tune_StageII(Palp,k2,m,museq(i),sampleX,y,bound,p,delta);
end
[~,s2ind] =sort(mybic2);
optmu = museq(s2ind(1));
%norm of penalized additive function estimation under optimal tuning parameter
[~,Pbeta,Pnbeta,~] = tune_StageII(Palp,k2,m,optmu,sampleX,y,bound,p,delta);
optind3 = find(Pnbeta<bound);
%normbeta = mean(Pbeta(:,optind3).^2);
%linearInd = find(abs(normbeta-1/3)<bound);
normbeta = mean(Pbeta(:,optind3));
linearInd = find(abs(normbeta-1/2)<bound); 
end

