function sol = fun_MOR_combined(par,x0,kx0,kz0)
%FUN_MOR_SHEAR Integrate shear part of growth rate along streamlines
%   Detailed explanation goes here

C= par.theta1 - sin(par.theta1)*cos(par.theta1);

s_react=par.n*par.beta_s;
s_shear=2*par.lambda_s/(4/3+par.zeta_r);

theta0=atan(x0);
y0 = [0;0;kx0;kz0];
options = odeset('AbsTol',1e-6,'RelTol',1e-9);
sol = ode45(@ode_theta,[theta0 par.theta1],y0,options);

    function dydt = ode_theta(t,y)
        theta=t;
        r = Theta(theta0) / ( cos(theta0) * Theta(theta) );
        x = r*sin(theta);
        z = -r*cos(theta);
        [~,~,gradu,~,theta_e] = fun_MOR_base_mantle(par,x,z);
        w0_rel = fun_MOR_base_magma(par,x,z+1e-14);
        kx=y(3);
        kz=y(4);
        k_theta=atan(kz./kx);
        
        dydt=zeros(size(y));
        dydt(1) = w0_rel*r/Theta(theta)*s_react*(cos(k_theta))^2;
        dydt(2) = s_shear*abs(sin(theta))/C/(Theta(theta))*cos(2*(k_theta-theta_e));
        dydt(3:4) = -r/Theta(theta)*[gradu{1,1}*kx + gradu{2,1}*kz;...
                                gradu{1,2}*kx + gradu{2,2}*kz];
    end

    function Theta = Theta(theta) 
        Theta = (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;
    end


end



