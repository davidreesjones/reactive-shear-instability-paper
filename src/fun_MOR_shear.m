function sol = fun_MOR_shear(par,x0)
%FUN_MOR_SHEAR Integrate shear part of growth rate along streamlines
%   Detailed explanation goes here

C= par.theta1 - sin(par.theta1)*cos(par.theta1);

s1=2*par.lambda_s/(4/3+par.zeta_r);

theta0=atan(x0);
s0 = 0;
options = odeset('AbsTol',1e-6,'RelTol',1e-9);
sol = ode45(@ode_theta,[theta0 par.theta1],s0,options);

    function dydt = ode_theta(t,y)
        theta=t;
        r = Theta(theta0) / ( cos(theta0) * Theta(theta) );       
        dydt = s1*sin(abs(theta))/C/(Theta(theta));
    end

    function Theta = Theta(theta) 
        Theta = (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;
    end


end



