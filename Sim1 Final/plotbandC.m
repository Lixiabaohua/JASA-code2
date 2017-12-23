function  plotbandC(Dat, ealp, ebet, sig2, res, kC, kA, m, optm1, optN, Q, delta )
% comprehensioned  confidence bands
%% data
x=Dat(:,6:7);
u=Dat(:,8);
T=size(Dat,1);
% estimated asymptotic variancee
[Valp, Vbet]=estVar(T,sig2, ealp, ebet, Dat(:,8), Dat(:,6:7), m, kC, kA, delta );
%%  bootstrap confidence bands
Balp0=zeros(T,Q);   Balp1=zeros(T,Q); Balp2=zeros(T,Q);
Bbeta1=zeros(T,Q);   Bbeta2=zeros(T,Q);
for i=1:Q
   y= ealp(:, 1)+ealp(:, 2) .* ebet(:, 1)+ealp(:, 3) .* ebet(:, 2)+res .* normrnd(0, 1, T, 1);
   [halp,hbeta] =   Spest( optN, kC, kA, optm1, m, x, u, y, delta);
   Balp0(:, i)   = halp(:, 1);
   Balp1(:, i)   = halp(:, 2);
   Balp2(:, i)   = halp(:, 3);
   Bbeta1(:, i) = hbeta(:,1);
   Bbeta2(:, i) = hbeta(:, 2);
end
%% 
[sx,ind]=sort(x);
[salp0,~] = sort(Balp0,2);
[salp1,~] = sort(Balp1,2) ;
[salp2,~] = sort(Balp2,2) ;
[sbeta1,~] = sort(Bbeta1,2);
[sbeta2,~] = sort(Bbeta2,2) ;
i2 =ceil(Q*0.025); i3 = ceil(Q*0.975);
Balp=[salp0(:,i2) salp0(:,i3) salp1(:,i2) salp1(:,i3) salp2(:,i2) salp2(:,i3)];
Bbet1  =  [ sbeta1(:,i2) sbeta1(:,i3)];
Bbet2  = [ sbeta2(:,i2) sbeta2(:,i3)];
%% asymptotic confidence bands based on CLT
%for varying-coefficient function
% pointwise confidence bands 
LAlp1=norminv(0.975,0,1)*(Valp.^(1/2))*(kC/T)^(1/2);
lowAlp1=ealp-LAlp1;
upAlp1=ealp+LAlp1;
% simultaneous confidence bands
I1=(2*log(kC+1))^(-1);
I2=log(log(kC+1))+log(4*pi);
d1=1-I1*(log(0.025/m)+0.5*I2);
LAlp2=(2*m*log(kC+1))^(1/2)*d1*(Valp.^(1/2))*(kC/T)^(1/2);
lowAlp2 = ealp-LAlp2;
upAlp2  = ealp+LAlp2;
% for additive functions
% pointwise confidence bands
LBet1=norminv(0.975,0,1)*(Vbet.^(1/2))*(kA/T)^(1/2);
lowBet1= ebet-LBet1;
upBet1 = ebet+LBet1;
% simultaneous confidence bands
l1=(2*log(kA+1))^(-1);
l2=log(log(kA+1))+log(4*pi);
d2=1-l1*(log(0.025/m)+0.5*l2);
LBet2=(2*m*log(kA+1))^(1/2)*d2*(Vbet.^(1/2))*(kA/T)^(1/2);
lowBet2=ebet-LBet2;
upBet2=ebet+LBet2;
%% plot pointwise and simultaneous based on CLT
plot(u, ealp(:,1), 'r:',   u, Dat(:,1), 'k-',   u, lowAlp1(:,1), 'b--',    u, upAlp1(:,1), 'b--', ...
                                                          u, lowAlp2(:,1), 'b-.',   u, upAlp2(:,1), 'b-.')
xlim([-0.05 1.05]) 
ylim([-2.5 3.5])
xlabel('Rescaled Time')
ylabel('\alpha_{0}')


plot(u, ealp(:,2), 'r:',   u, Dat(:,2), 'k-',   u, lowAlp1(:,2), 'b--',    u,upAlp1(:,2), 'b--', ...
                                                         u, lowAlp2(:,2), 'b-.',    u, upAlp2(:,2), 'b-.')
xlim([-0.05 1.05]) 
xlabel('Rescaled Time')
ylabel('\alpha_{1}')   

 plot(u, ealp(:,3), 'r:',   u, Dat(:,3), 'k-',   u, lowAlp1(:,3), 'b--',   u, upAlp1(:,3), 'b--',...
                                                          u, lowAlp2(:,3), 'b-.',  u, upAlp2(:,3), 'b-.')
xlim([-0.05 1.05]) 
ylim([-2 3.5])
xlabel('Rescaled Time')
ylabel('\alpha_{2}')

plot(sx(:,1), ebet(ind(:,1),1),       'r:',                 sx(:,1),    Dat(ind(:,1), 4),         'k-', ...
       sx(:,1), lowBet1(ind(:,1),1),  'b--',              sx(:,1),    upBet1(ind(:,1),1),     'b--',...    
       sx(:,1), lowBet2(ind(:,1),1),  'b-.',   sx(:,1),    upBet2(ind(:,1),1),     'b-.' )
xlim([-2.2 2.2])
xlabel('X_{1}')
ylabel('\beta_{1}')

 plot(sx(:,2),  ebet(ind(:,2),2),      'r:',                sx(:,2),   Dat(ind(:,2), 5),           'k-', ...
        sx(:,2),  lowBet1(ind(:,2),2), 'b--',              sx(:,2),   upBet1(ind(:,2),2),      'b--', ...
        sx(:,2),  lowBet2(ind(:,2),2), 'b-.',              sx(:,2),   upBet2(ind(:,2),2),       'b-.' ) 
xlim([-1.2 1.3])
xlabel('X_{2}')
ylabel('\beta_{2}')
%% pointwise confidence based on CLT and bootstrap method
plot(u, ealp(:,1), 'r:',   u, salp0(:, T/2),  'k-',  u,  lowAlp1(:,1), 'b--',              u, upAlp1(:,1),   'b--', ...
                                                                u,  Balp(:,1),      'magenta-.',   u, Balp(:,2),       'magenta-.' )
xlim([-0.05 1.05])   
ylim([-2 2.5])
xlabel('Rescaled Time')
ylabel('\alpha_{0}')

plot(u, ealp(:,2), 'r:',   u, salp1(:, T/2),  'k-',  u, lowAlp1(:,2), 'b--',              u, upAlp1(:,2),  'b--', ...
                                                                u, Balp(:,3),       'magenta-.',  u, Balp(:,4),       'magenta-.') 
 xlim([-0.05 1.05])   
 ylim([-0.2 2])    
xlabel('Rescaled Time')
ylabel('\alpha_{1}')

plot(u, ealp(:,3), 'r:',   u, salp2(:, T/2),  'k-',  u, lowAlp1(:,3), 'b--',             u, upAlp1(:,3),  'b--', ...
                                                                u, Balp(:,5),      'magenta-.',  u, Balp(:,6),       'magenta-.')
xlim([-0.05 1.05])   
xlabel('X_{1}')
ylabel('\beta_{1}')
 xlabel('Rescaled Time')
ylabel('\alpha_{2}')

plot(sx(:,1), ebet(ind(:,1),1),      'r:',                sx(:,1),   sbeta1(ind(:,1),T/2),    'k-',    ...  
       sx(:,1), lowBet1(ind(:,1),1), 'b--',             sx(:,1),  upBet1(ind(:,1),1),       'b--',  ...
       sx(:,1), Bbet1(ind(:,1),1),    'magenta-.',  sx(:,1),  Bbet1(ind(:,1),2),         'magenta-.')  
 xlim([-2.2 2.2])
xlabel('X_{1}')
ylabel('\beta_{1}')

plot(sx(:,2),  ebet(ind(:,2),2),      'r:',                  sx(:,2),  sbeta2(ind(:,2),T/2),     'k-', ...
       sx(:,2),  lowBet1(ind(:,2),2), 'b--',               sx(:,2),   upBet1(ind(:,2),2),       'b--', ...
       sx(:,2),  Bbet2(ind(:,2),1),     'magenta-.',   sx(:,2),   Bbet2(ind(:,2),2),         'magenta-.' )
ylim([-4.6 1.2])
xlim([-1.3 1.3])
xlabel('X_{2}')
ylabel('\beta_{2}')
end

