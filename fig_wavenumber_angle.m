%% Choose parameter values to plot
v_theta_e=[0 pi/8 pi/4 3*pi/8 pi/2];
%v_theta_e=[pi/4 3*pi/8 ];
vS=[0 0.5  1 2 10 ];
% vS=[0 0.25 0.5 0.75 1 2 10 100];

%% Set-up arrays and labels
NJ=length(v_theta_e);
NI=length(vS);
str_y=cell(1,NJ);
str_y{1}='$\theta_e=0$';
str_y{2}='$\theta_e=\pi/8$';
str_y{3}='$\theta_e=\pi/4$';
str_y{4}='$\theta_e=3\pi/8$';
str_y{5}='$\theta_e=\pi/2$';

[THETA,PSI]=meshgrid(linspace(-pi/2,pi/2),linspace(0, pi));
KX=cos(THETA).*cos(PSI);
KY=cos(THETA).*sin(PSI);
KZ=sin(THETA);

f1 = figure(2); clf; 
t1 = tiledlayout(NJ,NI);
load_colormap

%% Calculate and make plots
for J=1:NJ
    theta_e=v_theta_e(J);
    for I=1:NI
        S=vS(I);
        sig_t=(cos(THETA)).^2+S*( ((cos(THETA).*cos(PSI)).^2-(sin(THETA)).^2)*cos(2*theta_e) +cos(PSI).*sin(2*THETA)*sin(2*theta_e));
        nexttile 
        contourf(THETA,PSI,sig_t/max(abs(sig_t(:))),39,'LineStyle','none')
        colormap(cmap)
        caxis([-1 1])
        box off
        ax=gca;
        axis off
        axis equal
        hold on
        plot(pi/2,pi/2,'x m',LW{:},MS{:})
        xc=0.5*atan(2*S*sin(2*theta_e)/(1+2*S*cos(2*theta_e)));
        if S>-cos(2*theta_e)
            plot(0,pi/2,'x m',LW{:},MS{:})
            if xc>=0
                plot(xc,0,'s m',LW{:},MS{:})
                plot(xc-pi/2,0,'o m',LW{:},MS{:})
            else
                plot(xc+pi/2,0,'s m',LW{:},MS{:})
                plot(xc,0,'o m',LW{:},MS{:})
            end
        else
            plot(0,pi/2,'s m',LW{:},MS{:})
            if xc>=0
                plot(xc,0,'x m',LW{:},MS{:})
                plot(xc-pi/2,0,'o m',LW{:},MS{:})
            else
                plot(xc+pi/2,0,'x m',LW{:},MS{:})
                plot(xc,0,'o m',LW{:},MS{:})
            end
        end
        if J==1; title(['$S=',num2str(S),'$'],TX{:},FSl{:}); end        
        if I==1; ylabel(str_y{J},TX{:},FSl{:}); end
        ax.Title.Visible='on';
        ax.YLabel.Visible='on';
    end
end

%% Formatting and output
t1.TileSpacing='none';
t1.Padding='compact';

title_str=strcat('Increasing matrix shear $\rightarrow$');
title(t1,title_str,TX{:},FSl{:})
ylabel(t1,'Angle of maximum extension',TX{:},FSl{:})
cbar=colorbar;
cbar.Location='southoutside';
cbar.TickLabelInterpreter=TX{2};
drawnow
f1.WindowStyle='normal';
f1.Units='centimeters';
f1.Position=[0 0 20 20];
print(strcat(path_spec,'wavenumber_angle'),output_format)
