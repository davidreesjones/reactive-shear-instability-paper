load_par
C= par.theta1 - sin(par.theta1)*cos(par.theta1);
f1=figure(8); clf;
t1=tiledlayout(4,1);
[x,z]=meshgrid(linspace(0,1.5,401),linspace(-1,0,201));
Nx=100;
Nth=200;
vx0=linspace(x(1,2),tan(par.theta1) -0.001,Nx);
vtheta0=atan(vx0);
vkx0=[0,1/sqrt(2),1,1/sqrt(2)];
vkz0=[1,1/sqrt(2),0,-1/sqrt(2)];
for J=1:4
    
    
    Ms_kx=zeros(Nx+1,Nth);
    Ms_kz=zeros(Nx+1,Nth);
    Mx=zeros(Nx+1,Nth);
    Mz=zeros(Nx+1,Nth);
    Theta = @(theta) (theta.*cos(theta)-sin(theta)*(cos(par.theta1))^2)/C;
    for I=1:Nx
        %vtheta_i=linspace(vtheta0(I),pi/2-alpha_MOR,Nth);
        vtheta_i=logspace(log10(vtheta0(I)),log10(par.theta1),Nth);
        vtheta_i(end)=par.theta1;
        vtheta_i(1)=vtheta0(I);
        vr_i=Theta(vtheta0(I))/cos(vtheta0(I))./Theta(vtheta_i);
        Mx(I,:)=vr_i.*sin(vtheta_i);
        Mz(I,:)=-vr_i.*cos(vtheta_i);
        kx0=vkx0(J);
        kz0=vkz0(J);
        sol = fun_MOR_wavenumber(par,vx0(I),kx0,kz0);
        Ms_kx(I,:)=deval(sol,vtheta_i,1);
        Ms_kz(I,:)=deval(sol,vtheta_i,2);
    end
    I=Nx+1;
    z_i=linspace(-1,0,Nth);
    Mz(I,:)=z_i;
    Ms_kx(I,:)=kx0+zeros(size(z_i));
    Ms_kz(I,:)=kz0+zeros(size(z_i));
    
    sq_kx = griddata(Mx,Mz,Ms_kx,x,z);
    sq_kz = griddata(Mx,Mz,Ms_kz,x,z);
    sq_kx(z>-x*tan(pi/2-par.theta1))=nan;
    sq_kz(z>-x*tan(pi/2-par.theta1))=nan;
    sq_ktheta = atan2(sq_kz,sq_kx);
    sq_kabs = sqrt(sq_kx.^2+sq_kz.^2);
    
    sq_kx(z>-x*tan(pi/2-par.theta1))=nan;
    sq_kz(z>-x*tan(pi/2-par.theta1))=nan;
    
    nexttile
    h=contourf(x,z,sq_kabs,20,'LineStyle','None');
    caxis([0 2])
    
    colormap(cmap)
    hold on
        plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
    [Mxc,Mzc]=meshgrid([0.1:0.15:1.4],[-1:0.15:0.1]);
    len=0.06;
    for I=1:numel(Mxc)
        xc=Mxc(I);
        zc=Mzc(I);
        if zc<-xc*tan(pi/2-par.theta1)
            k_theta = interp2(x,z,sq_ktheta,xc,zc);
            [h_k,h_k_perp] = wavenumber_plot(gca,LW,par,xc,zc,len,k_theta);
        end
    end
    if J==1
        cbar=colorbar;
        cbar.Location= 'eastoutside';
        cbar.LineWidth=LW{2};
        cbar.FontSize=FSs{2};
        cbar.TickLabelInterpreter=TX{2};
        cbar.Limits=[0 2];
        cbar.Ticks=0:0.5:2;
        set(get(cbar,'Title'),'String','$k$',TX{:},FSs{:})
    end
    axis equal
    box on
    ax=gca;
    ax.XTick=[0:0.5:1.5];
    ax.YTick=[-1:0.5:0];
    ax.FontSize=FSs{2};
    ax.TickLabelInterpreter=TX{2};
    ax.LineWidth=LW{2};
    switch J
        case 1
            %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
            text(0.85,-0.1,'(a): $\theta_0=\pi/2$',TX{:},FSm{:});
        case 2
            %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
            text(0.85,-0.1,'(b): $\theta_0=\pi/4$',TX{:},FSm{:});
        case 3
            %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
            text(0.85,-0.1,'(c): $\theta_0=0$',TX{:},FSm{:});
        case 4
            %title('$s_\mathrm{shear}$',TX{:},FSm{:});
            text(0.85,-0.1,'(d): $\theta_0=-\pi/4$',TX{:},FSm{:});
            xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
            xlab.Position(2)=-1.02;
    end
    ylab=ylabel('$\frac{z}{H}$',TX{:},FSm{:},'Rotation',0);
    ylab.Position(1:2)=[-0.1 -0.3];
    t1.TileSpacing='none';
t1.Padding='compact';
drawnow

f1.WindowStyle='normal';
f1.Units=UN{2};
f1.Position=[0 0 9 20];
print(strcat(path_spec,'MOR_wavenumber'),output_format)
    
end