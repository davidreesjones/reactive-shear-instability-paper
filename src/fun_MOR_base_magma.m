function w0_rel = fun_MOR_base_magma(par,x,z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n=par.n;
C= par.theta1 - sin(par.theta1)*cos(par.theta1);

w0_rel = par.F_max^(1-1/n)...
            *( ((1+x.^2).^(-1) - (cos(par.theta1))^2).*(1+z) / C).^(1-1/n)...
            *par.Q_rel^(1/n);
        


end

