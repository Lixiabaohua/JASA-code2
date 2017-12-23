function Ini= StepIpre( k, m, sampleX, y,delta )
%segment length 100
T=length(y);         %sample size
%generate B-spline in Step I
b1 = rspline(sampleX(:,1),sampleX(:,1), m,k);   
b2 = rspline(sampleX(:,2),sampleX(:,2), m,k);   
%design matrix in Step I estiamtion
B = [ones(T,1) b1 b2 ];
%Step I estimation: segment 
s1 = zeros(size(b1,2),1);  s2 = zeros(size(b2,2),1);    
n1 = floor(T/100); I1 =100;
for l = 1 : n1
      low = (l - 1 ) * I1 +1 ;
      up = I1 * l ;
      bx = B(low : up, : ) ;
      by = y( low : up) ;
      tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
      s1 = s1 + tem3(2:(size(b1,2)+1));
      s2 = s2 + tem3((size(b1,2)+2):(2*size(b1,2)+1));
end  
low = n1*100 + 1;
up =T;
bx = B(low : up, : ) ;
by = y( low : up) ;
tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
s1 = s1 + tem3(2:size(b1,2)+1);
s2 = s2 + tem3(size(b1,2)+2:2*size(b1,2)+1);
%initial estimate for additive function
hbeta1 =  b1 * s1;
hbeta2 =  b2 * s2;
%centering
Inibeta1 = hbeta1 -mean(hbeta1);
Inibeta2 = hbeta2 -mean(hbeta2);
Ini = [Inibeta1 Inibeta2];
end

