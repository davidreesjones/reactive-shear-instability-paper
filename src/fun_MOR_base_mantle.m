function [ux,uz,gradu,dot_gamma,theta_e] = fun_MOR_base_mantle(par,x,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

C= par.theta1 - sin(par.theta1)*cos(par.theta1);

r=sqrt(x.^2+z.^2);
theta=atan(x./(-z));

ux=( x.*z./r.^2 - atan(x./z) ) / C;
uz=(z.^2./r.^2 - (cos(par.theta1))^2 ) / C;

gradu = cell(2,2);
gradu{1,1}=-2.*z.*x.^2/C./r.^4;
gradu{1,2}=2*x.^3/C./r.^4;
gradu{2,1}=-2*x.*z.^2/C./r.^4;
gradu{2,2}=+2.*z.*x.^2/C./r.^4;            

dot_gamma = abs(x)./r.^2/C;
theta_e=theta-pi/4*sign(theta); 


end

