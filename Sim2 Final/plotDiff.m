function  plotDiff(Dat,para,theta )
% plot regression function under null and alternative hypothesis

% fix rescaled time at u=0.5 and plot curves about x
T=size(Dat,1);
xq = quantile(Dat(:,4),[0.1,0.99]);
mygrid = xq(1):0.1:xq(2);
beta1 = 3*sin(pi*mygrid/2) - (1-mygrid).*mygrid-para(3);
term0 = repmat(Dat(T/2,1),length(mygrid),length(theta));
term1 = repmat(Dat(T/2,2),length(mygrid),length(theta));
term3 = repmat(theta,length(mygrid),1);
term4 = repmat(mygrid'+0.5,1,length(theta));
f=term0+term1.*repmat(beta1',1,length(theta))+term3.*term4;
f1 = f(:,1); f2 = f(:,2); f3 = f(:,3); f4=f(:,4); f5 = f(:,5); f6 = f(:,6);
  
plot(mygrid, f1, 'k-', mygrid , f2, 'b--', mygrid , f3, 'r--', mygrid , f4, 'g--', mygrid , f5, 'm--',...
       mygrid ,f6,'c--','LineWidth',0.8)
end

