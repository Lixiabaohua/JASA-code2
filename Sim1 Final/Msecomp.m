function [ MseU, MseO, MseA, MseV,Msey, optknot] = Msecomp( Q,Dat,repY,m,m1seq,Nseq,delta)
%generate Table for given data(sample size)
alpMat = Dat(:,1:3);
betaMat = Dat(:,4:5);
xMat = Dat(:,6:7);
t = Dat(:,8);
y = Dat(:,9);
avey = Dat(:,1) +Dat(:,2).*Dat(:,4)+Dat(:,3).*Dat(:,5);
%sample size and covariates dimension
T=size(xMat,1);  
%choose optimal knots for varying-coefficient additive model
kseq = ceil(0.5*T^(1/5)):ceil(2*T^(1/5));
[optknot,~] = myknot_vca(kseq,m,m1seq,Nseq,xMat,t,y,delta) ;
optkC=optknot(1); optkA = optknot(2); optm1= optknot(3); optN=optknot(4);

Umse = zeros(Q,size([alpMat betaMat],2));   Omse = zeros(Q,size([alpMat betaMat],2));
Amse = zeros(Q,size(alpMat,2));                   Vmse = zeros(Q,size(alpMat,2));  
Uy = zeros(T,Q) ;  Oy = zeros(T,Q) ; 
Vy=zeros(T,Q) ;  Ay = zeros(T,Q) ;
  for i = 1: Q
    %three-step spline estimation
    y = repY(:,i);
    [Ualp,Ubeta,~] =Spest( optN,optkC, optkA, optm1,m, xMat, t, y,delta);
    Uy(:,i) = Ualp(:,1)+Ualp(:,2).*Ubeta(:,1) +Ualp(:,3).*Ubeta(:,2);
    %oracle estimation
    [Oalp,Obeta] = Oracle_est(alpMat,betaMat,t,xMat,y,optkC,optkA,m,delta ) ;
    Oy(:,i)= Oalp(:,1)+Oalp(:,2).*Obeta(:,1) +Oalp(:,3).*Obeta(:,2);
    % misspecified additive model
    [Valp,Abeta] =Mis_est( t,xMat,y,optkC,optkA,m,delta ) ;
    Vy(:,i) = Valp(:,1)+Valp(:,2).*Dat(:,6)+Valp(:,3).*Dat(:,7);
    Ay(:,i) = Abeta(:,1)+Abeta(:,2)+Abeta(:,3);
    %mse for different estimators
    Umse(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp Ubeta]).^2));
    Omse(i,:) =  sqrt(mean(([alpMat betaMat]-[Oalp Obeta]).^2));
    Amse(i, :) = sqrt(mean((Abeta-[alpMat(:,1) betaMat]).^2));
    Vmse(i,:) = sqrt(mean((Valp - alpMat).^2));
  end
 %MseUy = sqrt(mean((Uy-repmat(avey,1,Q)).^2));
 %MseOy = sqrt(mean((Oy-repmat(avey,1,Q)).^2));
 %MseVy = sqrt(mean((Vy-repmat(avey,1,Q)).^2));
 %MseAy = sqrt(mean((Ay-repmat(avey,1,Q)).^2));
 %Msey = [mean(MseUy) std(MseUy) mean(MseOy) std(MseOy) ...
  %             mean(MseVy) std(MseVy) mean(MseAy) std(MseAy)];
 %MseU = [mean(Umse) mean(MseUy); std(Umse) std(MseUy)];
 %MseO = [mean(Omse) mean(MseOy); std(Omse) std(MseOy)];
 %MseV = [mean(Vmse) mean(MseVy); std(Vmse) std(MseVy)];
 %MseA = [mean(Amse) mean(MseAy);  std(Amse) std(MseAy)];
 %MseU = [mean(Umse) ;std(Umse) ];
 %MseO = [mean(Omse) ; std(Omse)];
 %MseV = [mean(Vmse) ; std(Vmse) ];
 %MseA = [mean(Amse);  std(Amse) ];
 
 
 MseUy = sqrt(mean((Uy-repmat(avey,1,Q)).^2));
 MseOy = sqrt(mean((Oy-repmat(avey,1,Q)).^2));
 MseVy = sqrt(mean((Vy-repmat(avey,1,Q)).^2));
 MseAy = sqrt(mean((Ay-repmat(avey,1,Q)).^2));
 Msey = [mean(MseUy) std(MseUy) mean(MseOy) std(MseOy) ...
               mean(MseVy) std(MseVy) mean(MseAy) std(MseAy)];
 MseU = [mean(Umse) mean(MseUy); std(Umse) std(MseUy)];
 MseO = [mean(Omse) mean(MseOy); std(Omse) std(MseOy)];
 MseV = [mean(Vmse) mean(MseVy); std(Vmse) std(MseVy)];
 MseA = [mean(Amse) mean(MseAy);  std(Amse) std(MseAy)];
end


