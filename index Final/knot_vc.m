function optk= knot_vc( kseq,m,Sx,u,Sy,delta )
% knots selection for varying-coefficient model
[T,p] = size(Sx); 
bicseq1 = zeros(length(kseq),1);
for i =1:length(kseq)
    [~,res] = vcm(Sx, u,Sy, m, kseq(i), delta ) ;
    rss = sum(res.^2);
    para = (p+1)*(m+kseq(i));
    bicseq1(i)= log(rss)+para*log(T)/T;
end
optk = find(bicseq1==min(bicseq1));
end