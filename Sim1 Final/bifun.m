function [ fun1,fun2] = bifun( stdpara,x1,x2,t)
%compute product of two univariate functions
%[ X1P, TP1 ] = meshgrid( x1, t ) ;
%[ X2P, TP2 ] =meshgrid( x2, t) ;
%alp1= (2*TP1.*sin(2*pi*TP1) +1)/stdpara(1);
%alp2 = (3*(1-TP2).*cos(2*pi*TP2)+1)/stdpara(2);
%beta1 = 3*sin(pi* X1P/2) - (1- X1P).*X1P - stdpara(3);
%beta2 = 3*cos(pi* X2P/2) +1.8*sin(pi* X2P/3)-stdpara(4);
%fun1 = alp1.*beta1;
%fun2 = alp2.*beta2;
%the following is equivalent to the above method
alp1 = (1.3*t.*sin(2*pi*t) +1)/stdpara(1);
alp2=(2*sin(1.5*pi*t)-1.2*(t-0.5).*(1-t)+1)/stdpara(2);

beta1 = 3*sin(pi* x1/2) - (1- x1).*x1 - stdpara(3);
beta2 = 2*cos(pi*x2/2) +1.8*sin(pi*x2/3)-stdpara(4);
fun1 = kron(beta1,alp1);
fun1 = reshape(fun1,length(t),length(x1));
fun2 = kron(beta2,alp2);
fun2 = reshape(fun2,length(t),length(x2));
end

