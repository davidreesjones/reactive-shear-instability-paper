function [h_ext,h_comp] = strain_plot(ax,plot_opts,par,xc,zc,len)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[~,~,~,~,theta_e] = fun_MOR_base_mantle(par,xc,zc);
dx=len/2*cos(theta_e);
dz=len/2*sin(theta_e);
h_comp=plot(ax,[xc-dz,xc+dz],[zc+dx,zc-dx],'Color','r',plot_opts{:});
h_ext=plot(ax,[xc-dx,xc+dx],[zc-dz,zc+dz],'Color','k',plot_opts{:});
end

