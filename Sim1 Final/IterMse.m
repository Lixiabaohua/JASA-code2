function [MseU,Msey,optknot ] =IterMse( Q,Dat,repY,m,m1seq,Nseq,delta )
% Mse comparison of multi step iteration:
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

Umse0 = zeros(Q,size([alpMat betaMat],2));  
Umse1 = zeros(Q,size([alpMat betaMat],2));  
Umse2 = zeros(Q,size([alpMat betaMat],2));  
Umse3 = zeros(Q,size([alpMat betaMat],2)); 
%Umse4 = zeros(Q,size([alpMat betaMat],2)); 


Uy0 = zeros(T,Q) ;  Uy1 = zeros(T,Q) ;  Uy2 = zeros(T,Q) ;
Uy3 = zeros(T,Q) ;  %Uy4 = zeros(T,Q) ;


  for i = 1: Q
    %three-step spline estimation
    y = repY(:,i);
    % three-step estimation
    [Ualp0,Ubeta0,~] =Spest( optN,optkC, optkA, optm1,m, xMat, t, y,delta);
    Uy0(:,i) = Ualp0(:,1)+Ualp0(:,2).*Ubeta0(:,1) +Ualp0(:,3).*Ubeta0(:,2);
    
    % one time iteration(2 steps)
    [Ualp1,Ubeta1] =IterSpest(xMat,t,y,Ubeta0,optkA,optkC,m,delta);
    Uy1(:,i)= Ualp1(:,1)+Ualp1(:,2).*Ubeta1(:,1) +Ualp1(:,3).*Ubeta1(:,2);
    
    % two times iteration(3 steps)
    [Ualp2,Ubeta2] =IterSpest(xMat,t,y,Ubeta1,optkA,optkC,m,delta);
    Uy2(:,i)= Ualp2(:,1)+Ualp2(:,2).*Ubeta2(:,1) +Ualp2(:,3).*Ubeta2(:,2);
    
    % three times iteration(4 steps)
    [Ualp3,Ubeta3] =IterSpest(xMat,t,y,Ubeta2,optkA,optkC,m,delta);
    Uy3(:,i)= Ualp3(:,1)+Ualp3(:,2).*Ubeta3(:,1) +Ualp3(:,3).*Ubeta3(:,2);
    
    % four times iteration(5 steps)
    %[Ualp4,Ubeta4] =IterSpest(xMat,t,y,Ubeta3,optkA,optkC,m,delta);
    %Uy4(:,i)= Ualp4(:,1)+Ualp4(:,2).*Ubeta4(:,1) +Ualp4(:,3).*Ubeta4(:,2);
    
    %mse for different estimators
    Umse0(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp0 Ubeta0]).^2));
    Umse1(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp1 Ubeta1]).^2));
    Umse2(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp2 Ubeta2]).^2));
    Umse3(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp3 Ubeta3]).^2));
    %Umse4(i,:)  = sqrt(mean(([alpMat betaMat]-[Ualp4 Ubeta4]).^2));
  end
  
 MseUy0 = sqrt(mean((Uy0-repmat(avey,1,Q)).^2));
 MseUy1 = sqrt(mean((Uy1-repmat(avey,1,Q)).^2));
 MseUy2 = sqrt(mean((Uy2-repmat(avey,1,Q)).^2));
 MseUy3 = sqrt(mean((Uy3-repmat(avey,1,Q)).^2));
 %MseUy4 = sqrt(mean((Uy4-repmat(avey,1,Q)).^2));
 
 Msey = [mean(MseUy0) std(MseUy0) mean(MseUy1)  std(MseUy1)...
     mean(MseUy2) std(MseUy2)  mean(MseUy3) std(MseUy3)];
               
 MseU = [(mean(Umse0))' (std(Umse0))' (mean(Umse1))' (std(Umse1))' (mean(Umse2))' (std(Umse2))'...
                (mean(Umse3))' (std(Umse3))' ];
 

end

