 function [optknot, optbic] = myknot(kseq,m,m0,I1,n1,I2,n2,sampleX,u,y,delta )
%choose optimal knots
[T,p] = size(sampleX); 
bicseq1 = zeros(length(kseq),length(kseq)); 
for i =1:length(kseq)
    for j = 1:length(kseq)
    inib = StepIest( kseq(j), m0, I1, n1,I2,n2,sampleX, y,delta ) ;
    [~,~,~,Ures,~] = Spest(inib, kseq(i), kseq(j), m, sampleX, u, y,delta) ;
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

