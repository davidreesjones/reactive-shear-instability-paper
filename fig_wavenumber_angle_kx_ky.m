%% Choose parameter values to plot
v_theta_e=[0 pi/8 pi/4 3*pi/8 pi/2];
%v_theta_e=[pi/4 3*pi/8 ];
vS=[0 0.4  1 2 10 ];
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

[THETA,PSI]=meshgrid(linspace(0,pi/2,40),linspace(0, 2*pi,40));
KX=cos(THETA).*cos(PSI);
KY=cos(THETA).*sin(PSI);
KZ=sin(THETA);

f1 = figure(4); clf; 
f1.WindowStyle='normal';
f1.Units='centimeters';
f1.Position=[0 0 26 23];
t1 = tiledlayout(NJ+1,NI+1);
load_colormap

%% Calculate and make plots
for J=1:NJ
    theta_e=v_theta_e(J);
    nexttile(1+(J-1)*(NI+1))
    box off
    axis on
    axis equal
    
    ax=gca;
    x = linspace(-1,1,101);
    y = linspace(-1,1,101);
    c = [0.15 0.3 0.45];
    c = [-c,0,c];
    [X Y] = meshgrid(x,y);
    Xp = X*cos(theta_e) + Y*sin(theta_e);
    Yp = Y*cos(theta_e) - X*sin(theta_e);
    psi = Xp.*Yp;
    h = contour(x,y,psi,c,'k','linewidth',1); hold on
    plot([-1 1],[0 0],'-k');
    plot([0 0],[-1 1],'-k');
    set(gca,'color','none','visible','off','xlim',[-1.5 1.5],'ylim',[-1.5 1.5]);
 
    for I=1:NI
        S=vS(I);
        sig_t=(cos(THETA)).^2+S*( ((cos(THETA).*cos(PSI)).^2-(sin(THETA)).^2)*cos(2*theta_e) +cos(PSI).*sin(2*THETA)*sin(2*theta_e));
        nexttile(1+I+(J-1)*(NI+1)) 
        contourf(KX,KY,sig_t/max(abs(sig_t(:))),39,'LineStyle','none')
        colormap(cmap)
        caxis([-1 1])
        box off
        ax=gca;
        axis off
        axis equal
        hold on
        %plot(0,0,'x m',LW{:},MSl{:})
        xc=0.5*atan(2*S*sin(2*theta_e)/(1+2*S*cos(2*theta_e)));
        if S~=0
        if S>-cos(2*theta_e)
            %plot(0,1,'x m',LW{:},MSl{:})
            if xc>=0
                plot(cos(xc),0,'s m',LW{:},MSl{:})
                % plot(cos(xc-pi/2),0,'o m',LW{:},MSl{:})
            else
                plot(cos(xc+pi/2),0,'s m',LW{:},MSl{:})
                % plot(cos(xc),0,'o m',LW{:},MSl{:})
            end
        else
            plot(0,1,'s m',LW{:},MSl{:})
            if xc>=0
                % plot(cos(xc),0,'x m',LW{:},MSl{:})
                %plot(cos(xc-pi/2),0,'o m',LW{:},MSl{:})
            else
                % plot(cos(xc+pi/2),0,'x m',LW{:},MSl{:})
                %plot(cos(xc),0,'o m',LW{:},MSl{:})
            end
        end
        end
        if J==1; title(['$S=',num2str(S),'$'],TX{:},FSl{:}); end        
        if I==1; ylabel(str_y{J},TX{:},FSl{:}); end
        ax.Title.Visible='on';
        ax.YLabel.Visible='on';
    end
end


%% plot with labels
nexttile((NI+1)*NJ+2)
[MKX,MKY]=meshgrid(linspace(-1,1,40),linspace(-1,1,40));
h= contourf(MKX,MKY,MKX.^2+MKY.^2,[1 1],'k',LW{:});
ax=gca;
axis off
axis equal
hold on
box off


nexttile((NI+1)*NJ+3)
h= contourf(KX,KY,sig_t/max(abs(sig_t(:))),39,'LineStyle','none','Visible','off');

cbar=colorbar;
cbar.Location='north';
cbar.TickLabelInterpreter=TX{2};
cbar.Label.String='Relative growth rate';
cbar.Label.FontSize=FSl{2};
cbar.FontSize=FSl{2};
cbar.Label.Interpreter=TX{2};
colormap(cmap)
caxis([-1 1])
axis off
axis equal
hold on
        
%% Formatting and output
t1.TileSpacing='none';
t1.Padding='compact';

title_str=strcat('Increasing matrix shear $\rightarrow$');
title(t1,title_str,TX{:},FSl{:})
yl=ylabel(t1,'$\leftarrow$ Angle of maximum extension',TX{:},FSl{:});

drawnow
for J=1:NJ
    theta_e=v_theta_e(J);
    nexttile(1+(J-1)*(NI+1))
    arrow_plot(gca,[0 0] ,[0.5*cos(theta_e),0.5*sin(theta_e)]);
    arrow_plot(gca,[0 0] ,[-0.5*cos(theta_e),-0.5*sin(theta_e)]);
    arrow_plot(gca,[-0.8*sin(theta_e),0.8*cos(theta_e)],[-0.3*sin(theta_e),+0.3*cos(theta_e)]);
    arrow_plot(gca,[0.8*sin(theta_e),-0.8*cos(theta_e)],[0.3*sin(theta_e),-0.3*cos(theta_e)]);
    text(1.1,0,'$x$',TX{:},FSm{:})
    text(0,1.25,'$z$',TX{:},FSm{:})

end
nexttile((NI+1)*NJ+2)
ar=arrow_plot(gca,[0 0] ,[0,1]); ar.LineStyle='-'; ar.LineWidth=LW{2};
ar=arrow_plot(gca,[0 0] ,[1,0]); ar.LineStyle='-'; ar.LineWidth=LW{2};
text(0.45,0.2,'$k_x$',TX{:},FSl{:})
text(0.05,0.7,'$k_y$',TX{:},FSl{:})      
print(strcat(path_spec,'wavenumber_angle_kx_ky'),output_format)
