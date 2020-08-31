f1=figure(16); clf
t1=tiledlayout(3,3);
ax=cell(3,3);
[x,z]=meshgrid(linspace(0.001,1.5,101),linspace(-0.99,-0.001,81));
options = optimset('Display','off','TolX',1e-6);

for J=1:3
    switch J
        case 1
            load_par
            par.Q_rel=par.Q_rel*10;
        case 2
            load_par
        case 3
            load_par
            par.Q_rel=0.4e4;
    end
    
    w0_rel = fun_MOR_base_magma(par,x,z);
    [ux,uz,gradu,dot_gamma,theta_e] = fun_MOR_base_mantle(par,x,z);
    S=2*par.lambda_s/(par.n*par.beta_s*(4/3+par.zeta_r))*dot_gamma./w0_rel;
    S(z>-x*tan(pi/2-par.theta1))=nan;
    theta_max=0.5*atan2(2*S.*sin(2*theta_e),1+2*S.*cos(2*theta_e));
    
    
    
    M_k_theta_c=zeros(size(z))+nan;
    M_k_theta_end=zeros(size(z))+nan;
    M_s_total_end=zeros(size(z))+nan;
    for I=1:numel(x)
        xc=x(I); zc=z(I);
        if zc<-xc*tan(pi/2-par.theta1)
            k_theta_c=fminbnd(@(k_theta) -fun_MOR_s_total_interior(par,xc,zc,k_theta),-pi/2,pi/2,options);
            [s_total,sol] = fun_MOR_s_total_interior(par,xc,zc,k_theta_c);
            M_k_theta_c(I)=k_theta_c;
            kz=sol.y(4,end);
            if abs(imag(kz))<1e-14 %floating point issues
                kz=real(kz);
            end
            kx=sol.y(3,end);
            if abs(imag(kx))<1e-14 %floating point issues
                kx=real(kx);
            end
            M_k_theta_end(I)=atan2(kz,kx);
            M_s_total_end(I)=s_total;
        end
    end
    
    nexttile
    ax{J,1}=gca; hold on;
    contourf(x,z,M_s_total_end,200,'LineStyle','None')
    caxis([-40 40])
    colormap(cmap)
    if J==3
        cbar=colorbar;
        cbar.Limits(1)=0;

        cbar.Location= 'southoutside';
        cbar.LineWidth=LW{2};
        cbar.FontSize=FSs{2};
        cbar.TickLabelInterpreter=TX{2};
    end
    
    nexttile
    ax{J,2}=gca; hold on;
    contourf(x,z,M_k_theta_end,200,'LineStyle','None')
    caxis([-pi/4 pi/4])
    colormap(cmap)
    %title('$\theta_\mathrm{max}$',TX{:},FSm{:})
    if J==3
        cbar=colorbar;
        cbar.Location= 'southoutside';
        cbar.LineWidth=LW{2};
        cbar.FontSize=FSs{2};
        cbar.TickLabelInterpreter=TX{2};
    end
    
    nexttile
    ax{J,3}=gca; hold on;
    h=contourf(x,z,(theta_max),200,'LineStyle','None');
    caxis([-pi/4 pi/4])
    colormap(cmap)
    if J==3
        cbar=colorbar;
        cbar.Location= 'southoutside';
        cbar.LineWidth=LW{2};
        cbar.FontSize=FSs{2};
        cbar.TickLabelInterpreter=TX{2};
    end
    h=contour(ax{J,2},x,z,M_s_total_end,[7 7],'-- m',LW{:});
end

for J=1:9
    axis(ax{J},'equal')
    box(ax{J}, 'on');
    hold on;
    plot(ax{J},[0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
    ax{J}.XTick=[0:0.5:1.5];
    ax{J}.YTick=[-1:0.5:0];
    ax{J}.FontSize=FSs{2};
    ax{J}.TickLabelInterpreter=TX{2};
    ax{J}.LineWidth=LW{2};
    set(ax{J} ,'Layer', 'Top')
end
title(ax{1,1},'Amplitude: $s_\mathrm{total}$',TX{:},FSm{:})
title(ax{1,2},'Global:  $\theta$',TX{:},FSm{:})
title(ax{1,3},'Local: $\theta_\mathrm{max}$',TX{:},FSm{:})
txt1=text(ax{1,1},0.4,-0.1,'(a): fast melt velocity',TX{:},FSm{:});
txt2=text(ax{2,1},0.4,-0.1,'(b): ref. melt velocity',TX{:},FSm{:});
txt3=text(ax{3,1},0.4,-0.1,'(c): slow melt velocity',TX{:},FSm{:});

t1.TileSpacing='none';
t1.Padding='compact';
drawnow

f1.WindowStyle='normal';
f1.Units=UN{2};
f1.Position=[0 0 20 15];
print(strcat(path_spec,'MOR_theta_max'),output_format)