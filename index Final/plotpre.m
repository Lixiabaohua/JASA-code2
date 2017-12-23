function plotpre(y,pre,err, T)
% predict at the last 50 points
 
 u=2969:T;
 
 plot(u, err) 
 xlim([2965 3020])
 boxplot(err)
 

 plot(u,y(u),'k-o',u,pre,'b-.x')
 xlim([2965 3020])
 
end
