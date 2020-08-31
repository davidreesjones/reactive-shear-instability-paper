function [h_k,h_k_perp] = wavenumber_plot(ax,plot_opts,par,xc,zc,len,k_theta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dx=len/2*cos(k_theta);
dz=len/2*sin(k_theta);
h_k=plot(ax,[xc,xc+2*dx],[zc,zc+2*dz],'Color','k',plot_opts{:});
h_k_perp=plot(ax,[xc-dz,xc+dz],[zc+dx,zc-dx],'Color','r',plot_opts{:});
end

