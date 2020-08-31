load_par

C= par.theta1 - sin(par.theta1)*cos(par.theta1);
f1=figure(13); clf;
t1=tiledlayout(3,3);
[x,z]=meshgrid(linspace(0,1.5,401),linspace(-1,0,201));
Nx=100;
Nth=200;
vx0=linspace(x(1,2),tan(par.theta1) -0.001,Nx);
vtheta0=atan(vx0);
vkx0=0.963842158559942;
vkz0=sqrt(1-vkx0^2);
xf0=0.2;
thetaf0=atan(xf0);
vtheta=linspace(thetaf0,par.theta1);
r0=sec(thetaf0);
vr=r0*Theta(thetaf0)./Theta(vtheta);
vx=vr.*sin(vtheta);
vz=-vr.*cos(vtheta);

for J=1:3
    if J==1
        load_par
        par.zeta_r=10;
    elseif J==2
        load_par
        par.Q_rel=par.Q_rel*10;
    else
        load_par
        par.Q_rel=0.4e4;
    end
    
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
%         kx0=vkx0(J);
%         kz0=vkz0(J);
        options = optimset('Display','off','TolX',1e-6);
        theta_c=fminbnd(@(theta) -fun_MOR_s_total_end(par,vx0(I),cos(theta),sin(theta)),-pi/2,pi/2,options);
        sol = fun_MOR_combined(par,vx0(I),cos(theta_c),sin(theta_c));
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
        if J==1
            caxis([-20 20])
        elseif J==2
            caxis([-60 60])
        else
            caxis([-10 10])
        end
        colormap(cmap)
        hold on
        plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
        %plot(vx,vz,'g',LW{:})
        [Mxc,Mzc]=meshgrid([0.1:0.2:1.4],[-1:0.15:0.1]);
        len=0.05;
        for I=1:numel(Mxc)
            xc=Mxc(I);
            zc=Mzc(I);
            if zc<-xc*tan(pi/2-par.theta1)
                k_theta = interp2(x,z,sq_ktheta,xc,zc);
                [h_k,h_k_perp] = wavenumber_plot(gca,LW,par,xc,zc,len,k_theta);
            end
        end
        if K==1&&J==1                   
            text(0.4,-0.1,'(a): bulk viscosity high',TX{:},FSm{:});
        end
        if K==1&&J==2
            text(0.4,-0.1,'(b): melt velocity high',TX{:},FSm{:});
        end
        if K==1&&J==3
            text(0.4,-0.1,'(c): melt velocity low',TX{:},FSm{:});
        end
        if K==3
                cbar=colorbar;
                cbar.Location= 'eastoutside';
                cbar.LineWidth=LW{2};
                cbar.FontSize=FSs{2};
                cbar.TickLabelInterpreter=TX{2};
                if J==1
                    cbar.Limits=[-20 20];
                    cbar.Ticks=-20:10: 20;
                elseif J==2
                    cbar.Limits=[-60 60];
                    cbar.Ticks=-60:30: 60;
                else
                    cbar.Limits=[-10 10];
                    cbar.Ticks=-10:5: 10;
                end
        end
        axis equal
        box on
        ax=gca;
        ax.XTick=[0:0.5:1.5];
        ax.YTick=[-1:0.5:0];
        ax.FontSize=FSs{2};
        ax.TickLabelInterpreter=TX{2};
        ax.LineWidth=LW{2};
        if K==1
            ylab=ylabel('$\frac{z}{H}$',TX{:},FSm{:},'Rotation',0);
            ylab.Position(1:2)=[-0.1 -0.3];
        end
        xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
        xlab.Position(2)=-1.02;
    end
    title1.FontSize=12;
    title2.FontSize=12;
    title3.FontSize=12;
    t1.TileSpacing='compact';
    t1.Padding='compact';
    drawnow
    
    f1.WindowStyle='normal';
    f1.Units=UN{2};
    f1.Position=[0 0 20 14];
    print(strcat(path_spec,'MOR_growth_favoured_sensitivity'),output_format)
    
end