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

[THETA,PSI]=meshgrid(linspace(-pi/2,pi/2,20),linspace(0, pi,10));
KX=cos(THETA).*cos(PSI);
KY=cos(THETA).*sin(PSI);
KZ=sin(THETA);

f1 = figure(2); clf; 
f1.WindowStyle='normal';
f1.Units='centimeters';
f1.Position=[0 0 25 23];
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


%% plot with labels
nexttile((NI+1)*NJ+2)
h=plot([-pi/2 pi/2],[0 pi],'LineWidth',1,'Color','w');
box on
ax=gca;
axis on
axis equal
hold on
set(ax,'TickLabelInterpreter','latex')
set(ax,'XLim',[-pi/2 pi/2],'YLim',[0 pi])
set(ax,'XTick',[-pi/2 pi/2],'XTickLabel',{'$-\pi/2$','$\pi/2$'})
set(ax,'YTick',[0 pi],'YTickLabel',{'$0$','$\pi$'})
ylabel('$\psi$','Rotation',0,TX{:},FSl{:})
xlabel('$\theta$','Rotation',0,TX{:},FSl{:},'Position',[0 -0.4])
set(ax,'TickLabelInterpreter','latex',FS{:})
box on

nexttile((NI+1)*NJ+3)
h= contourf(THETA,PSI,sig_t/max(abs(sig_t(:))),39,'LineStyle','none','Visible','off');

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
      
print(strcat(path_spec,'wavenumber_angle'),output_format)
