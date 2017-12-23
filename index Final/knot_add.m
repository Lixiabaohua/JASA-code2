function [optknot,optbic]= knot_add(kseq,sampleX,y,m,delta )
%optimal knots for additive model
[T,p] = size(sampleX);
bic1 = zeros(1,length(kseq));
for i = 1:length(kseq)
    [~,rss,~]  = add_est(sampleX,y,m,kseq(i),delta);
    bic1(i) = log(rss)+(1+p*(kseq(i)+m))*log(T)/T;
end
[~,ind]=sort(bic1);
optknot = kseq(ind(1));
optbic = min(bic1);
end

