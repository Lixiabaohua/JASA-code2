function [ fun ] = deriveSCAD( x,lambda )
%compute the first derivative of SCAD penalization function 
fun = lambda*(x<=lambda)+max(3.7*lambda-x,0)*(x>lambda)/2.7;
end

