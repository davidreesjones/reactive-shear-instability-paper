load_par
f1=figure(7); clf;
t1=tiledlayout(2,2);
[x,z]=meshgrid(linspace(0,1.5,401),linspace(-1,0,201));
[ux,uz,gradu,dot_gamma,theta_e] = fun_MOR_base_mantle(par,x,z);
dot_gamma(z>-x*tan(pi/2-par.theta1))=nan;

nexttile
levels=-2.5:0.1:2.5;
h=contourf(x,z,log10(dot_gamma),levels,'LineStyle','None');
colormap(gca,cmap_BG)
cbar=colorbar;
cbar.Location= 'eastoutside';
cbar.LineWidth=LW{2};
cbar.FontSize=FSs{2};
cbar.TickLabelInterpreter=TX{2};
cbar.Limits=[-0.5 1];
cbar.Ticks=[-0.5 0 0.5 1 2];
for I=1:numel(cbar.TickLabels)
    cbar.TickLabels{I}=strcat('$10^{',num2str(str2num(cbar.TickLabels{I})),'}$');
end
caxis([-0.5 1])
ax=gca;
axis equal
box on;
hold on;
plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
[Mxc,Mzc]=meshgrid([0.1:0.1:1.4],[-0.9:0.1:0.1]);
len=0.06;
for I=1:numel(Mxc)
    xc=Mxc(I);
    zc=Mzc(I);
    if zc<-xc*tan(pi/2-par.theta1)
        strain_plot(ax,LW,par,xc,zc,len);
    end
end
ax.XTick=[0:0.5:1.5];
ax.YTick=[-1:0.5:0];
ax.FontSize=FSs{2};

%xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
%xlab.Position(2)=-1.02;
ylab=ylabel('$\frac{z}{H}$',TX{:},FSm{:},'Rotation',0);
ylab.Position(1:2)=[-0.1 -0.3];

title('Normalized strain rate $\dot{\gamma}_0/(U_0/H)$',TX{:},FSm{:})
text(1.39,-0.06,'(a)',TX{:},FSm{:});
ax.TickLabelInterpreter=TX{2};
ax.LineWidth=LW{2};
w0_rel = fun_MOR_base_magma(par,x,z);
[ux,uz,gradu,dot_gamma,theta_e] = fun_MOR_base_mantle(par,x,z);
w0_rel(z>-x*tan(pi/2-par.theta1))=nan;
dot_gamma(z>-x*tan(pi/2-par.theta1))=nan;
S=2*par.lambda_s/(par.n*par.beta_s*(4/3+par.zeta_r))*dot_gamma./w0_rel;
for J=1:3
    nexttile
    
    
    switch J
        case 1
            levels=0:5:130;
            h=contourf(x,z,(w0_rel),levels,'LineStyle','None');
            colormap(gca,cmap_Br)
        case 2
            levels=[0:10:250]*1e-5;
            h=contourf(x,z,(w0_rel/par.Q_rel).^(1/(par.n-1)),levels,'LineStyle','None');
            colormap(gca,cmap_Br)
        case 3
            h=contourf(x,z,log10(S),40,'LineStyle','None');
            hold on
            for x0=[0.01,0.1:0.1:1.3]
            h_stream = streamline_plot(gca,LW,'m',par,x0);
            if x0==0.01||x0==0.2
                h_stream.LineStyle='--';
            end
            end
            colormap(gca,cmap_BrBG)
            caxis([-1.5 1.5])
    end
    cbar=colorbar;
    cbar.Location= 'eastoutside';
    cbar.LineWidth=LW{2};
    cbar.FontSize=FSs{2};
    cbar.TickLabelInterpreter=TX{2};
    cbar.Limits=[levels(1) levels(end-2)];
    cbar.Ticks=levels(1:4:end);
    if J==3
        cbar.Limits=[-1 1];
        cbar.Ticks=-1:1:1;
        for I=1:numel(cbar.TickLabels)
            cbar.TickLabels{I}=strcat('$10^{',num2str(str2double(cbar.TickLabels{I})),'}$');
        end
        caxis([-1 1])
    end
    ax=gca;
    axis equal
    box on;
    hold on;
    plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
    ax.XTick=[0:0.5:1.5];
    ax.YTick=[-1:0.5:0];
    ax.FontSize=FSs{2};
    
    switch J
        case 1
            title('Normalized magma velocity: $w_0/U_0$',TX{:},FSm{:});
            text(1.39,-0.06,'(b)',TX{:},FSm{:});
        case 2
            title('Porosity: $\phi_0$',TX{:},FSm{:});
            text(1.39,-0.06,'(c)',TX{:},FSm{:});
            xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
            xlab.Position(2)=-1.02;
            ylab=ylabel('$\frac{z}{H}$',TX{:},FSm{:},'Rotation',0);
            ylab.Position(1:2)=[-0.1 -0.3];
        case 3
            title('Ratio of shear to reactive growth rate: $S$',TX{:},FSm{:});
            text(1.39,-0.06,'(d)',TX{:},FSm{:});
            xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
            xlab.Position(2)=-1.02;
    end
    ax.TickLabelInterpreter=TX{2};
    ax.LineWidth=LW{2};
end
t1.TileSpacing='compact';
t1.Padding='compact';
drawnow

f1.WindowStyle='normal';
f1.Units=UN{2};
f1.Position=[0 0 20 13];
print(strcat(path_spec,'MOR_base_state'),output_format)