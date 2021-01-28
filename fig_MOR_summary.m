options = optimset('Display','off','TolX',1e-6);
load_par
depth_scale=60; %km

vx0=[0.01,0.2,0.5];
ind=1:numel(vx0);
c_sol=cell(1,numel(vx0));
ax=cell(1,numel(vx0));
title_str=cell(1,numel(vx0));
f1=figure(15);
t1=tiledlayout(1,numel(vx0)+1);

for x0=vx0
    switch ind(x0==vx0)
        case 1
            prefix='near';
            lab='(a) $\,$';
        case numel(vx0)
            prefix='far';
            lab='(c) $\,$ '; %broken if not three panels
        otherwise
            prefix='mid';
            lab='(b) $\,$';
    end
    
    nexttile
    theta_c=fminbnd(@(theta) -fun_MOR_s_total_end(par,x0,cos(theta),sin(theta)),-pi/2,pi/2,options);
    sol = fun_MOR_combined(par,x0,cos(theta_c),sin(theta_c));
    vartheta=sol.x;
    theta0=atan(x0);
    r=Theta(theta0)/cos(theta0)./Theta(vartheta);
    z=-r.*cos(vartheta);
    x=r.*sin(vartheta);
    xf=x(end);
    title_str{x0==vx0}=strcat(prefix,': $x_1=', num2str(xf*depth_scale,'%.0f'),'$ km');
    s_r=sol.y(1,:);
    s_s=sol.y(2,:);
    Y=zeros(2,numel(vartheta));
    Y(1,:)=s_r;
    Y(2,:)=s_s;
    A=area(z*depth_scale,Y');
    A(1).FaceColor=map_BrBG(2,:);
    A(2).FaceColor=map_BrBG(end-1,:);
    
    view(90,270);
    title1=title(strcat(lab,title_str(x0==vx0)),TX{:},FSm{:});
    
    set(gca,'XLim',[-depth_scale,0],'YLim',[0 40])
    c_sol{x0==vx0}=sol;
    ax{x0==vx0}=gca;
    ylabel('$s$',TX{:},FSm{:});
    if x0==vx0(1)
        A(1).DisplayName='reaction';
        A(2).DisplayName='shear';
        l=legend;
        l.Location='southeast';
    end
    
end
nexttile
hold on;
colororder(newcolours([3,4,5]))

for x0=vx0
    sol=c_sol{x0==vx0};
    k_theta=atan2(sol.y(4,:),sol.y(3,:));
    vartheta=sol.x;
    theta0=atan(x0);
    r=Theta(theta0)/cos(theta0)./Theta(vartheta);
    z=-r.*cos(vartheta);
    plot(k_theta*180/pi,z*depth_scale,LW{:},'DisplayName',title_str{x0==vx0})
    set(gca,'YLim',[-depth_scale,0])
    ax{numel(vx0)+1}=gca;
    xlabel('$\theta$ (degrees)',TX{:},FSl{:});
    l2=legend;
    l2.Location='southwest';
        title1=title('(d) wavevector orientation',TX{:},FSl{:});

    
end
for K=1:numel(vx0)+1
    ax{K}.FontSize=FSm{2};
    ax{K}.TickLabelInterpreter=TX{2};
    ax{K}.LineWidth=LW{2};
    box(ax{K},'on')
end
xlabel(ax{1},'$\mathrm{depth (km)}$',TX{:},FSm{:});

l.FontSize=FSms{2};
l.Interpreter=TX{2};
l.ItemTokenSize(1)=10;

l2.FontSize=FSms{2};
l2.Interpreter=TX{2};
l2.ItemTokenSize(1)=10;

t1.TileSpacing='compact';
t1.Padding='compact';
drawnow

f1.WindowStyle='normal';
f1.Units=UN{2};
f1.Position=[0 0 20 6.5];
print(strcat(path_spec,'MOR_summary'),output_format)