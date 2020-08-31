load_par
C= par.theta1 - sin(par.theta1)*cos(par.theta1);
f1=figure(9); clf;
t1=tiledlayout(2,1);
[x,z]=meshgrid(linspace(0,1.5,401),linspace(-1,0,201));
[ux,uz,gradu,dot_gamma,theta_e] = fun_MOR_base_mantle(par,x,z);
dot_gamma(z>-x*tan(pi/2-par.theta1))=nan;

k_theta=0;
Nx=100;
Nth=200;
vx0=linspace(x(1,2),tan(par.theta1) -0.001,Nx);
vtheta0=atan(vx0);
Ms_shear=zeros(Nx+1,Nth);
Ms_react=zeros(Nx+1,Nth);
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
        
    sol = fun_MOR_shear(par,vx0(I));
    Ms_shear(I,:)=deval(sol,vtheta_i);
    
    sol = fun_MOR_react(par,vx0(I));
    Ms_react(I,:)=deval(sol,vtheta_i);
end

I=Nx+1;
z_i=linspace(-1,0,Nth);
Mz(I,:)=z_i;
n=par.n;
Ms_shear(I,:)=zeros(size(z_i));
s1=par.n*par.beta_s;
Ms_react(I,:)=s1*par.F_max^(1-1/n)*(C*(sin(par.theta1))^-2)^(1/n)*par.Q_rel^(1/n)*(1+z_i).^(2-1/n)*n/(2*n-1);

sq_shear = griddata(Mx,Mz,Ms_shear,x,z);
sq_react = griddata(Mx,Mz,Ms_react,x,z);


for J=1:2
    nexttile 
    ax=gca;
    axis equal
    box on;
    hold on;
    switch J
        case 1
            levels=0:1:25;
            h=contourf(x,z,sq_react,levels,'LineStyle','None');
            colormap(gca,cmap_Br)
            h2=contour(x,z,sq_react,[7 7],'-- m',LW{:})
            caxis([0 25])
        case 2
            levels=0:1:25;
            h=contourf(x,z,sq_shear,levels,'LineStyle','None');
            colormap(gca,cmap_BG)
            h2=contour(x,z,sq_shear,[7 7],'-- m',LW{:})
            caxis([0 25])
    end
    cbar=colorbar;
    cbar.Location= 'eastoutside';
    cbar.LineWidth=LW{2};
    cbar.FontSize=FSs{2};
    cbar.TickLabelInterpreter=TX{2};
    cbar.Limits=[levels(1) levels(end)];
    cbar.Ticks=levels(1:5:end);

    
    plot([0 1/tan(pi/2-par.theta1)],[0 -1],'k',LW{:});
    ax.XTick=[0:0.5:1.5];
    ax.YTick=[-1:0.5:0];
    ax.FontSize=FSs{2};
    
    switch J
        case 1
            title('$s_\mathrm{reaction}$',TX{:},FSm{:});
            text(1.39,-0.06,'(a)',TX{:},FSm{:});
        case 2
            title('$s_\mathrm{shear}$',TX{:},FSm{:});
            text(1.39,-0.06,'(b)',TX{:},FSm{:});
            xlab=xlabel('$\frac{x}{H}$',TX{:},FSm{:});
            xlab.Position(2)=-1.02;
    end
    ylab=ylabel('$\frac{z}{H}$',TX{:},FSm{:},'Rotation',0);
    ylab.Position(1:2)=[-0.1 -0.3];
    
    ax.TickLabelInterpreter=TX{2};
    ax.LineWidth=LW{2};
    set(ax,'Layer','top');
end
t1.TileSpacing='compact';
t1.Padding='compact';
drawnow

f1.WindowStyle='normal';
f1.Units=UN{2};
f1.Position=[0 0 10 13];
print(strcat(path_spec,'MOR_growth_upper'),output_format)