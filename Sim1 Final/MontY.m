function RepY = MontY( MyDat ,Q)
% generate Monte Carlo relication
% MyDat give the varying-coefficient component functions (Column 1:3) and additive
% component functions(Column 4:5)
T=size(MyDat,1);
RepY = zeros(T,Q);
t=MyDat(:,8);

for i = 1: Q
    RepY(:,i) =MyDat(:,1)+MyDat(:,2) .* MyDat(:,4) + MyDat(:,3) .*  MyDat(:,5)  + ...
          (0.7+(2-t)./(2+t)).*normrnd(0,0.4,T,1);
end
end

