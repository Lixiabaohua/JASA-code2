function [Valp, Vbeta]=estVar(T,sig,eAlp,eBeta,u,x,m,kC,kA,delta )
%construct the asymptotic variance of CLT
% sig2: estimated variance
% eAlp and eBeta are  estimated three-step spline estimation of component
% function
% u and x are sample od rescaled time and covariates

% asymptotic variance  for varying-coefficient functions
Bt= rspline(u, u ,m, kC); 

hatD1 = (kron(eBeta(:,1), ones(1,size(Bt,2)))).*Bt;
hatD2=(kron(eBeta(:,2), ones(1,size(Bt,2)))).*Bt;
hatD=[Bt hatD1 hatD2];

hatE1=hatD'*hatD/T;
temp1=pinv(hatE1+delta*eye(size(hatE1,2)));

Temp1=hatD.*repmat(sig.^(1/2),1,size(hatD,2));
Sig1=Temp1'*Temp1/T;

blkD1=kron(eye(T),ones(1,3*size(Bt,2)));

A0=[Bt zeros(T,2*size(Bt,2))];
A1=[zeros(T,size(Bt,2)) Bt zeros(T,size(Bt,2))];
A2=[zeros(T,2*size(Bt,2)) Bt];

LA0=repmat(A0,1,T).*blkD1;
LA1=repmat(A1,1,T).*blkD1;
LA2=repmat(A2,1,T).*blkD1;

prod1=kron(eye(T),temp1*Sig1*temp1);

RA0=reshape(A0',T*size(A0,2),1);
RA1=reshape(A1',T*size(A1,2),1);
RA2=reshape(A2',T*size(A2,2),1);

% vector of asymptotic variance
VarAlp0=LA0*prod1*RA0;
VarAlp1=LA1*prod1*RA1;
VarAlp2=LA2*prod1*RA2;


% asymptotic variance  for additive functions
b1 = rspline(x(:,1),x(:,1), m,kA);   
b2 = rspline(x(:,2),x(:,2), m,kA); 

hatZ1 =  (kron(eAlp(:,1), ones(1,size(b1,2)))).*b1;
hatZ2 =  (kron(eAlp(:,2), ones(1,size(b2,2)))).*b2;
hatZ  = [hatZ1 hatZ2 ];

hatE2=hatZ'*hatZ/T;
temp2=pinv(hatE2+delta*eye(size(hatE2,2)));

Temp2=hatZ.*repmat(sig.^(1/2),1,size(hatZ,2));
Sig2=Temp2'*Temp2/T;

blkD2=kron(eye(T),ones(1,2*size(b1,2)));

B1=[b1 zeros(T,size(b1,2))];
B2=[zeros(T,size(b1,2)) b2];

LB1=repmat(B1,1,T).*blkD2;
LB2=repmat(B2,1,T).*blkD2;

prod2=kron(eye(T),temp2*Sig2*temp2);

RB1=reshape(B1',T*size(B1,2),1);
RB2=reshape(B2',T*size(B2,2),1);


% vector of asymptotic variance
VarBet1=LB1*prod2*RB1;
VarBet2=LB2*prod2*RB2;

Vbeta=[VarBet1 VarBet2];
Valp=[VarAlp0 VarAlp1 VarAlp2 ];

end

