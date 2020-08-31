function [h_stream] = streamline_plot(ax,plot_opts,col,par,x0)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    C= par.theta1 - sin(par.theta1)*cos(par.theta1);
Theta = @(theta) (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;

    theta0=atan(x0);
    vtheta=logspace(log10(theta0),log10(par.theta1),200);
    vr = (Theta(theta0)/cos(theta0))./Theta(vtheta);
    vx = vr.*sin(vtheta);
    vz = -vr.*cos(vtheta);
    h_stream=plot(ax,vx,vz,'Color',col,plot_opts{:});
end

