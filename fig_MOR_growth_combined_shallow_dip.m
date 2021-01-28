load_par_shallow_dip
f1=figure(14); clf;
t1=tiledlayout(4,3);
[x,z]=meshgrid(linspace(0,3.5,401),linspace(-1,0,201));
Nx=100;
Nth=200;
vx0=linspace(x(1,2),tan(par.theta1) -0.001,Nx);
vtheta0=atan(vx0);
vkx0=[0,1/sqrt(2),1,1/sqrt(2)];
vkz0=[1,1/sqrt(2),0,-1/sqrt(2)];
for J=1:4
    
    
    Ms_kx=zeros(Nx+1,Nth);
    Ms_kz=zeros(Nx+1,Nth);
    Ms_s_react=zeros(Nx+1,Nth);
    Ms_s_shear=zeros(Nx+1,Nth);
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
        sol = fun_MOR_combined(par,vx0(I),kx0,kz0);
        Ms_s_react(I,:)=deval(sol,vtheta_i,1);
        Ms_s_shear(I,:)=deval(sol,vtheta_i,2);
        Ms_kx(I,:)=deval(sol,vtheta_i,3);
        Ms_kz(I,:)=deval(sol,vtheta_i,4);
    end
    I=Nx+1;
    z_i=linspace(-1,0,Nth);
    Mz(I,:)=z_i;
    Ms_kx(I,:)=kx0+zeros(size(z_i));
    Ms_kz(I,:)=kz0+zeros(size(z_i));
    
    sq_kx = griddata(Mx,Mz,Ms_kx,x,z);
    sq_kz = griddata(Mx,Mz,Ms_kz,x,z);
    sq_s_react = griddata(Mx,Mz,Ms_s_react,x,z);
    sq_s_shear = griddata(Mx,Mz,Ms_s_shear,x,z);
    sq_kx(z>-x*tan(pi/2-par.theta1))=nan;
    sq_kz(z>-x*tan(pi/2-par.theta1))=nan;
    sq_s_react(z>-x*tan(pi/2-par.theta1))=nan;
    sq_s_shear(z>-x*tan(pi/2-par.theta1))=nan;
    sq_ktheta = atan2(sq_kz,sq_kx);
    sq_kabs = sqrt(sq_kx.^2+sq_kz.^2);
    
    for K=1:3
        
        nexttile
        switch K
            case 1
                h=contourf(x,z,sq_s_react+sq_s_shear,20,'LineStyle','None');
                if J==1; title1=title('$s_\mathrm{total}$',TX{:},FS{:}); end
            case 2
                h=contourf(x,z,sq_s_react,20,'LineStyle','None');
                if J==1; title2=title('$s_\mathrm{reaction}$',TX{:},FS{:}); end
            case 3
                h=contourf(x,z,sq_s_shear,20,'LineStyle','None');
                if J==1; title3=title('$s_\mathrm{shear}$',TX{:},FS{:}); end
        end
        caxis([-20 20])
        
        colormap(cmap)
        hold on
        plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
        [Mxc,Mzc]=meshgrid([0.1:0.4:3.4],[-0.9:0.2:0.1]);
        len=0.1;
        for I=1:numel(Mxc)
            xc=Mxc(I);
            zc=Mzc(I);
            if zc<-xc*tan(pi/2-par.theta1)
                k_theta = interp2(x,z,sq_ktheta,xc,zc);
                [h_k,h_k_perp] = wavenumber_plot(gca,LW,par,xc,zc,len,k_theta);
            end
        end
        %if J==1
        %         cbar=colorbar;
        %         cbar.Location= 'east';
        %         cbar.LineWidth=LW{2};
        %         cbar.FontSize=FSs{2};
        %         cbar.TickLabelInterpreter=TX{2};
        %cbar.Limits=[0 2];
        %cbar.Ticks=0:0.5:2;
        %end
        axis equal
        box on
        ax=gca;
        ax.XTick=[0:3.5];
        ax.YTick=[-1:0];
        ax.FontSize=FSm{2};
        ax.TickLabelInterpreter=TX{2};
        ax.LineWidth=LW{2};
        if K==1
            switch J
                case 1
                    %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
                    text(2.,-0.2,'(a): $\theta_0=\pi/2$',TX{:},FSm{:});
                case 2
                    %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
                    text(2.0,-0.2,'(b): $\theta_0=\pi/4$',TX{:},FSm{:});
                case 3
                    %title('$s_\mathrm{reaction}$',TX{:},FSm{:});
                    text(2.0,-0.2,'(c): $\theta_0=0$',TX{:},FSm{:});
                case 4
                    %title('$s_\mathrm{shear}$',TX{:},FSm{:});
                    text(2.0,-0.2,'(d): $\theta_0=-\pi/4$',TX{:},FSm{:});
            end
            ylab=ylabel('$\frac{z}{H}$',TX{:},FSl{:},'Rotation',0);
            ylab.Position(1:2)=[-0.2 -0.65];
        end
        if J==4
            xlab=xlabel('$\frac{x}{H}$',TX{:},FSl{:});
            xlab.Position(1:2)=[1.5,-1.07];
            if K==2
                cbar=colorbar;
                cbar.Location= 'southoutside';
                cbar.LineWidth=LW{2};
                cbar.FontSize=FSm{2};
                cbar.TickLabelInterpreter=TX{2};
                cbar.Limits=[-20 20];
                cbar.Ticks=-20:5:20;
            end
        end
    end
    title1.FontSize=14;
    title2.FontSize=14;
    title3.FontSize=14;
    t1.TileSpacing='compact';
    t1.Padding='none';
    drawnow
    
    f1.WindowStyle='normal';
    f1.Units=UN{2};
    f1.Position=[0 0 20 11];
    print(strcat(path_spec,'MOR_growth_combined_shallow_dip'),output_format)
    
end