function [ optknot,optbic ] = Bivknot( kseq, m, sampleX, u, y, delta)
% knots selections of bivariate function estimations
% Dat: DGP
% m: order of B-spline
T = size(sampleX,1); 
bicseq1 = zeros(length(kseq),length(kseq)); 
for i =1:length(kseq)
    for j = 1:length(kseq)
    [~, res] = BivEst( sampleX, u, y, m, kseq(i), kseq(j), delta);
    rss = sum(res.^2);
    para = m + kseq(i) + 2 * (m + kseq(i)) * (m + kseq(j));
    bicseq1(i, j)= log(rss)+para*log(T)/T;
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

