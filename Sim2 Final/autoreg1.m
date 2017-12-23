function [ ts ] = autoreg1( n,rho,sigma  )
%Generate tvNAR(1) 2*n data, and discard the first n
%initial value ts0 = 0
ts = zeros( 2 * n , 1 ) ;
tim = linspace( 0, 1,  2 * n + 1 ) ;
t = tim( 2 : end) ;
ts( 1 ) = sigma* normrnd( 0, 1 ) ;
for i = 2 : ( 2 * n )
      ts( i ) = rho * t ( i  ) * ts( i - 1) +  sigma * normrnd( 0, 1 ) ;
end
ts( 1 : n ) =  [ ] ;
end

