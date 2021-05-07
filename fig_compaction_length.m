%% Choose parameters to plot
vS = [ 0.2 1]; % series of S
vtheta_e=[pi/8 3*pi/8];
vd = logspace(-2,2); % range of (k\delta)^2 to plot

%% Set-up
f=figure(5); clf;
load_colormap;
t=tiledlayout(2,1);
for theta_e=vtheta_e
    
    h=cell(numel(vS),3);
    nexttile
    hold on; box on;
    ax=gca;
    
    %% Calculate and plot
    for S = vS
        theta=0.5*atan2(2*S*sin(2*theta_e),(1+2*S*cos(2*theta_e)));
        sig_react=(1+vd*(cos(theta))^2)./(1+vd);
        sig_shear=S*cos(2*(theta-theta_e))*vd./(1+vd);
        h{S==vS,1}=plot(vd,sig_react);
        h{S==vS,2}=plot(vd,sig_shear);
        h{S==vS,3}=plot(vd,sig_react+sig_shear);
        set(ax,'XScale','log')
    end
    
    h{1,1}(1).DisplayName='Reactive';
    h{1,2}(1).DisplayName='Shear';
    h{1,3}(1).DisplayName='Total';
    
    for I=1:numel(vS)
        h{I,1}(1).Color='b';
        h{I,2}(1).Color='r';
        h{I,3}(1).Color='k';
        for J=1:3
            h{I,J}(1).LineWidth=LW{2};
            if I~=1
                h{I,J}(1).LineStyle='--';
                h{I,J}(1).HandleVisibility='off';
            end
        end
    end
    ylabel('$\mathrm{real}(\sigma)/\sigma_\mathrm{reaction} $',TX{:},FS{:})
    if theta_e==vtheta_e(end)
        xlabel('$(k\delta)^2$',TX{:},FS{:})
        set(ax,'YLim',[-0.5,1.5])
    end
    ax.TickLabelInterpreter=TX{2};
    set(ax,FSm{:})
    set(ax,'LineWidth',1)
    set(ax,'Layer','top')
    
end
l=legend;
l.Location="west";
l.Interpreter=TX{2};
l.ItemTokenSize(1)=20;
l_title=title(l,{'$S=0.2$ {(solid)}', '$S=1$ (dash)'});
l_title.Interpreter=TX{2};
l.FontSize=FSms{2};




% title(t,TX{:},FSs{:})
f.Units='centimeters';
f.Position=[0 0 10 13];

t.TileSpacing = 'compact';
t.Padding='compact';

print(strcat(path_spec,'compaction-length'),output_format)
