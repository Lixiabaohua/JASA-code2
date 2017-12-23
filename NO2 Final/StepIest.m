function inibeta= StepIest( k, m, I1, sampleX, y,delta )
% construct the initial estimation in Step I estimation (additive function)
%  k: interior knots number in initial estimation
% m: order of B-spline in initial estimation
% N: groups numbers in step I estimation
 T=length(y);         %sample size
%generate B-spline in Step I
b1 = rspline(sampleX(:,1),sampleX(:,1), m,k);   
b2 = rspline(sampleX(:,2),sampleX(:,2), m,k);   
b3=  rspline(sampleX(:,3),sampleX(:,3), m,k);  
%design matrix in Step I estiamtion
B = [ones(T,1) b1 b2 b3 ];
%Step I estimation: segment 
s1 = zeros(size(b1,2),1);  s2 = zeros(size(b2,2),1);    s3 = zeros(size(b3,2),1);  
N = T/I1;  
for l = 1 : N
      low = (l - 1 ) * I1 +1 ;
      up = I1 * l ;
      bx = B(low : up, : ) ;
      by = y( low : up) ;
      tem3 = pinv( bx' * bx +delta *eye(size(bx,2)))* bx'  * by ;
      s1 = s1 + tem3(2:(size(b1,2)+1));
      s2 = s2 + tem3((size(b1,2)+2):(2*size(b1,2)+1));
      s3 = s3 + tem3((2*size(b1,2)+2):(3*size(b1,2)+1));
end        
%initial estimate for additive function
hbeta1 =  b1 * s1;
hbeta2 =  b2 * s2;
hbeta3 =  b3 * s3;
%centering
Inibeta1 = hbeta1 -mean(hbeta1);
Inibeta2 = hbeta2 -mean(hbeta2);
Inibeta3 = hbeta3 -mean(hbeta3);
inibeta = [Inibeta1 Inibeta2 Inibeta3];
end

