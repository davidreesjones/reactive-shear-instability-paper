function [s_total,sol] = fun_MOR_s_total_interior(par,xc,zc,k_theta)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
C= par.theta1 - sin(par.theta1)*cos(par.theta1);
Theta = @(theta) (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;
rc=sqrt(xc^2+zc^2);
thetac=atan2(xc,-zc);
theta0=fzero(@(theta0) rc*Theta(thetac)-Theta(theta0)/cos(theta0),[0 par.theta1]);
r0=1/cos(theta0);
x0=r0*sin(theta0);
kx0=cos(k_theta);
kz0=sin(k_theta);
sol = fun_MOR_combined_interior(par,x0,kx0,kz0,thetac);
s_react=sol.y(1,end);
s_shear=sol.y(2,end);        
s_total = s_react+s_shear;

end

