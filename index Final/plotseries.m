function plotseries(Dat)
y = Dat(:,2);
x = Dat(:,1)/100;

plot(y);
title('High-low Range')
xlabel('(a)')
xlim([-100 3100])
set(gca,'XTickLabel',{'2000','2002','2004','2006','2008','2010','2012'})


plot(x);
title('Return')
xlim([-100 3100])
xlabel('(b)')
set(gca,'XTickLabel',{'2000','2002','2004','2006','2008','2010','2012'})


end

