function optknot=myknot_vca( kseq,m,m0seq,Nseq,sampleX,u,y,delta )
%choose optimal knots, groups and order of B-spline in intial estimation
%Nseq is optional of segment length
[knot1, bic1] = myknot(kseq,m,m0seq(1),Nseq(1),sampleX,u,y,delta );
[knot2, bic2] = myknot(kseq,m,m0seq(1),Nseq(2),sampleX,u,y,delta );
[knot3, bic3] = myknot(kseq,m,m0seq(1),Nseq(3),sampleX,u,y,delta );
[knot4, bic4] = myknot(kseq,m,m0seq(2),Nseq(1),sampleX,u,y,delta );
[knot5, bic5] = myknot(kseq,m,m0seq(2),Nseq(2),sampleX,u,y,delta );
[knot6, bic6] = myknot(kseq,m,m0seq(2),Nseq(3),sampleX,u,y,delta );
[~,ind1] = sort([bic1 bic2 bic3 bic4 bic5 bic6]);
expr = ind1(1);
switch expr
    case 1
         optknot = [knot1 m0seq(1) Nseq(1)];
    case 2
        optknot =  [knot2 m0seq(1) Nseq(2)];
    case 3
        optknot = [knot3 m0seq(1) Nseq(3)];
    case 4
        optknot = [knot4 m0seq(2) Nseq(1)];
    case 5
        optknot = [knot5 m0seq(2) Nseq(2)];  
    case 6
        optknot = [knot6 m0seq(2) Nseq(3)];  
end
end

