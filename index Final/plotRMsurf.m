function plotRMsurf(Inib, xMat, y, t, kC, kA, m, delta)
%plot 3-dimensional surface plot
%% generate grid points
tseq = linspace( quantile(t, 0.02), quantile( t, 0.98), 20);
x1seq = linspace( quantile( xMat(:,1), 0.02), quantile( xMat(:,1), 0.98 ), 30) ;
x2seq = linspace( quantile( xMat(:,2), 0.03), quantile( xMat(:,2), 0.97 ), 30) ;
[xm1, tm1]=meshgrid(x1seq, tseq);
[ ~, ~, coefC, coefA, a, c] = RM( Inib, kC, kA, m, xMat, t, y, delta ) ;
[alp, beta] = pred_RMest(x1seq', x2seq', tseq', xMat, t, coefC, coefA, a, c, m, kC, kA ) ;
%% compute the values of surfaces
efun1 = kron(beta(:, 1), alp(:, 2));
efun1 = reshape(efun1,length(tseq),length(x1seq));
%% plot surfaces
mesh(xm1, tm1, efun1 ) ;
ylabel('u') ;
xlabel('Y_{t-1}') ;
zlabel('m_{1}') ;
box on 
end


