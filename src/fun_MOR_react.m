function sol = fun_MOR_react(par,x0)
%FUN_MOR_REACT Integrate reactive part of growth rate along streamlines
%   Detailed explanation goes here

C= par.theta1 - sin(par.theta1)*cos(par.theta1);

s1=par.n*par.beta_s;

theta0=atan(x0);
s0 = 0;
options = odeset('AbsTol',1e-6,'RelTol',1e-9);
sol = ode45(@ode_theta,[theta0 par.theta1],s0,options);

    function dydt = ode_theta(t,y)
        theta=t;
        r = Theta(theta0) / ( cos(theta0) * Theta(theta) );
        x = r * sin(theta);
        z = -r * cos(theta);
        w0_rel = fun_MOR_base_magma(par,x,z+1e-14);
        if imag(w0_rel)~=0
            x
            z
            w0_rel
            pause
        end
        dydt = w0_rel*r/Theta(theta)*s1;
    end

    function Theta = Theta(theta)
        Theta = (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;
    end


end



