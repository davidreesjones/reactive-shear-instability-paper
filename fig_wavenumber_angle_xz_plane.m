%% Choose range of theta_e to plot
vtheta_e=pi/2*[0.25 0.75]; 

%% Set-up
f=figure(3); clf;
load_colormap
theta=linspace(-pi,pi,1e3);
t=tiledlayout(numel(vtheta_e),1);
ax=cell(numel(vtheta_e),1);

%% Calculate and plot
for theta_e = vtheta_e
    nexttile
    ax{theta_e == vtheta_e}=gca;
    hold on; box on;
    A1 = area([theta_e - 3*pi/4, theta_e - pi/4],[1 1],-1);
    A2 = area([theta_e + pi/4, theta_e + 3*pi/4],[1 1],-1);
    vS = [ 0 0.2 1 10 100];
    str=num2str(8*theta_e/pi);
    str=str(~strcmp(num2str(8*theta_e/pi),'1'));
    h=plot([theta_e theta_e],[-1 1],'--k',LW{:},'DisplayName',strcat('$\theta_e=',str,'\pi/8$'),HV{:});
    set( ax{theta_e == vtheta_e},'ColorOrderIndex', 1);
    for S = vS
        sig_t=(cos(theta)).^2+S*cos(2*(theta-theta_e));
        h=plot(theta,sig_t/max(sig_t));
        h.DisplayName=strcat('$S=',num2str(S),'$');
        h.LineWidth=1;
        xc=0.5*atan(2*S*sin(2*theta_e)/(1+2*S*cos(2*theta_e)));
        if S>-cos(2*theta_e)
            if xc>=0
                hmark_max=plot(xc,1,'s m',MS{:});
                hmark_min=plot(xc-pi/2,((cos(xc-pi/2)).^2+S*cos(2*(xc-pi/2-theta_e)))/max(sig_t),'o m',MS{:});
            else
                hmark_max=plot(xc+pi/2,1,'s m',MS{:});
                hmark_min=plot(xc,((cos(xc)).^2+S*cos(2*(xc-theta_e)))/max(sig_t),'o m',MS{:});
            end
            hmark_min.LineWidth=1;
            hmark_min.HandleVisibility='off';
            hmark_max.LineWidth=1;
            hmark_max.HandleVisibility='off';
        else
            if xc>=0
                hmark_saddle=plot(xc,1,'x m',MS{:});
                hmark_min=plot(xc-pi/2,((cos(xc-pi/2)).^2+S*cos(2*(xc-pi/2-theta_e)))/max(sig_t),'o m',MS{:});
            else
                hmark_saddle=plot(xc+pi/2,1,'x m',MS{:});
                hmark_min=plot(xc,((cos(xc)).^2+S*cos(2*(xc-theta_e)))/max(sig_t),'o m',MS{:});
            end
            hmark_saddle.LineWidth=1;
            hmark_saddle.HandleVisibility='off';
            hmark_min.LineWidth=1;
            hmark_min.HandleVisibility='off';
        end
      
        
    end
    plot([-pi pi],[0 0],'-k',LW{:},HV{:})
    if theta_e == vtheta_e(1)
        l=legend;
        l.Location="northwest";
        l.Interpreter='latex';
        l.ItemTokenSize(1)=20;
        l.FontSize=FSs{2};
    end
    lab=text(theta_e+0.05,-1+0.1,'$\theta_e$');
    set(lab,FSm{:},TX{:})
    A1.LineStyle='none';
    A1.FaceColor=[0 0 1];
    A1.FaceAlpha=0.1;
    A1.HandleVisibility='off';
    A2.LineStyle='none';
    A2.FaceColor=[0 0 1];
    A2.FaceAlpha=0.1;
    A2.HandleVisibility='off';
    
end
    xlabel('$\theta$',TX{:},FS{:})

t.TileSpacing = 'compact';
t.Padding='compact';
for I=1:numel(ax)
    ax{I}.TickLabelInterpreter=TX{2};
    set(ax{I},'XLim',[-pi/2 pi/2],'YLim',[- 1 1])
    set(ax{I},FSm{:})
    set(ax{I},LW{:})
    set(ax{I},'XTick',[-pi/2, -3*pi/8 -pi/4, -pi/8, 0, pi/8 pi/4, 3*pi/8, pi/2],'YTick',[- 1 -0.5 0 0.5 1]);
    set(ax{I},'Layer','top')
    set(ax{I},'XTickLabel',{'$-\pi/2$', '$-3\pi/8$','$-\pi/4$','$-\pi/8$','$0$','$\pi/8$','$\pi/4$', '$3\pi/8$', '$\pi/2$'});
end
%title(t,'Normalized growth rate for modes restricted to $x$--$z$ plane',FSm{:},TX{:})
ylabel(t,'Normalized growth rate: $\mathrm{real}(\sigma) / \max\limits_\theta \mathrm{real}(\sigma)$ ',FSm{:},TX{:})

f.Units='centimeters';
f.Position=[0 0 10 11];

print(strcat(path_spec,'wavenumber_angle_xz_plane'),output_format)