function s_total = fun_MOR_s_total_end(par,vx0,kx0,kz0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sol = fun_MOR_combined(par,vx0,kx0,kz0);
s_react=deval(sol,par.theta1,1);
s_shear=deval(sol,par.theta1,2);        
s_total = s_react+s_shear;
end

