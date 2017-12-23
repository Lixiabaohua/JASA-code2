function out = mse_Pre(Q,Dat,resY, tn, kC,kA,m1,m,gnum,delta )
%prediction RMSE of component functions and response
%based on Q Monte Carlo replications
% Dat: sample 
% resY: replicated response
% tn:  dropped sample number to build up model estimation(length of
% predection)
% gnum: segment length
% m1,m: order of B-spline matrix
% kC,kA: interior knots number
T=size(Dat,1);    
ryn = zeros(Q,1); ryMVn = zeros(Q,1);  ryMAn = zeros(Q,1);
for i =1:Q 
y = resY(:,i);
%three-step estimation
ryn(i)= Sp_test(T,Dat(:,1:8), y, tn, kC,kA,m1,m,gnum,delta ) ;
%Misspecified estimation
[ryMVn(i),ryMAn(i)] = mis_test( T, Dat(:,1:8), y, tn, kC,kA,m,delta) ;
end
eyn  = [mean(ryn) std(ryn)];
eMVyn = [mean(ryMVn) std(ryMVn)] ;
eMAyn = [mean(ryMAn) std(ryMAn)] ;
out = [eyn eMVyn eMAyn];
end


