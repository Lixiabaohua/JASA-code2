function [optknot, optbic] = myknot(kseq,m,m1,ncase,sampleX,u,y,delta )
% choose optimal knots
% kseq: potential knot
% m1: order of B-spline in Step I estimation(fixed)
% m: order of B-spline in Step II and Step III estimation(fixed)
% ncase: segment length
[T,p] = size(sampleX); 
bicseq1 = zeros(length(kseq),length(kseq)); 
for i =1:length(kseq)
    for j = 1:length(kseq)
    [~,~,Ures] = Spest(ncase,kseq(i), kseq(j), m1, m, sampleX, u, y,delta) ;
    rss = sum(Ures.^2);
    para = (p+1)*(m+kseq(i))+p*(m+kseq(j));
    bicseq1(i,j)= log(rss)+para*log(T)/T;
    end
end
optind1 = find(bicseq1 == min(min(bicseq1)));
colknot1 = ceil(optind1/length(kseq));
rowknot1 = mod(optind1,length(kseq));
if rowknot1 ==0
    rowknot1 =length(kseq);
end
optknot = [kseq(rowknot1) kseq(colknot1)];
optbic = min(min(bicseq1));
end

