function [ myAdd,myVary,myAll] = Mytest( Aterms,Vterms,NormA,NormV,bound )
% identify additive terms, varying-coefficient terms and overall model
% Aterms: terms of additive term
% Vterms: terms of varying-coefficient term
% NormA: norm of additive component functions
% NormV: norm of varying-coefficient functions
  myaddT=0;   myaddU=0;    myaddO=0;
  myvaryT =0; myvaryU=0;    myvaryO=0; 
  myallT=0;     myallU=0;       myallO=0;
% aditive term identification
    if NormV(3)< bound && length(Aterms) == 1
       myaddT = 1 ; 
    end
   if NormV(3)< bound && length(Aterms) >=2
       myaddU =1 ; 
   end
   if   NormV(3)>= bound 
       myaddO=  1 ; 
   end
  %varying-coefficient term identification
   if  NormA(4)<bound && length(Vterms) == 1 
      myvaryT =   1 ;   
   end
   if    NormA(4)<bound && length(Vterms) >=2 
       myvaryU =  1 ;
   end
  if  NormA(4)>= bound
     myvaryO =   1 ; 
  end
 %overall model identification
   if  length(Aterms)==1 && length(Vterms)==1 && NormA(4)<bound && NormV(3)< bound
       myallT=1;
   elseif NormA(4)<bound && NormV(3)< bound && max(length(Aterms),length(Vterms))>=2
       myallU = 1 ;
   else
      myallO= 1;
   end
   myAdd=[myaddT myaddU myaddO];
   myVary=[myvaryT myvaryU myvaryO];
   myAll=   [myallT  myallU  myallO];
end

