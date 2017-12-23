function plotpre(u,fit,pre,err)
%compare  of one-step ahead prediction
 plot(u, err)
 
 boxplot(err)
 
 plot(u,fit(u),'k-o', u, pre,'b-.x')
 
end

