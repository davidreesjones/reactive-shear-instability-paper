load_par
C= par.theta1 - sin(par.theta1)*cos(par.theta1);
f1=figure(11); clf;
t1=tiledlayout(2,2);

vktheta0=linspace(-pi/2,pi/2,200);
vkx0=cos(vktheta0);
vkz0=sin(vktheta0);

% vkx0=[0,1/sqrt(2),1,1/sqrt(2)];
% vkz0=[1,1/sqrt(2),0,-1/sqrt(2)];
% vktheta0=atan(vkz0./vkx0);
vsol=cell(1,numel(vkx0));
vs_react=zeros(1,numel(vkx0));
vs_shear=zeros(1,numel(vkx0));
vkx=zeros(1,numel(vkx0));
vkz=zeros(1,numel(vkx0));
vx0=[0.01, 0.2];
ax=cell(2,2);
for K=1:2
    
    for J=1:numel(vkx0)
        kx0=vkx0(J);
        kz0=vkz0(J);
        sol = fun_MOR_combined(par,vx0(K),kx0,kz0);
        vsol{J}=sol;
        vs_react(J)=deval(sol,par.theta1,1);
        vs_shear(J)=deval(sol,par.theta1,2);
        vkx(J)=deval(sol,par.theta1,3);
        vkz(J)=deval(sol,par.theta1,4);
    end
    vktheta=atan(vkz./vkx);
    for J=numel(vktheta0)-1:-1:1
        if abs(vktheta(J)-vktheta(J+1))>pi/4            
            vktheta(J)=vktheta(J)-pi; %avoids jumps associated with nonuniquness of atan
        end
    end
    [~,I]=max(vs_react+vs_shear);
    vs_total=vs_react+vs_shear;
    I_0=find(vs_total>=vs_total(I)*0.8); %find region within 80% of max
    [~,I_L]=min(vktheta0(I_0));
    [~,I_R]=max(vktheta0(I_0));
    I_L=I_0(I_L); % fix indexing
    I_R=I_0(I_R);
    vktheta_L=vktheta(I_L);
    vktheta_R=vktheta(I_R);
    A_X=[-pi/2 -pi/2 vktheta0(I_L) vktheta0(I_L) vktheta0(I_R) vktheta0(I_R)];
    A_Y=[vktheta_R vktheta_L vktheta_L -pi -pi vktheta_R];
    for K2=1:2
        nexttile
        ax{K,K2}=gca;
        hold on
        box on
        switch K2
            case 1
                A=area(vktheta0(I_0),vs_total(I_0),-100);
            case 2
                A=fill(A_X,A_Y,'c');
        end
        plot([-pi/2 pi/2],[0 0],'Color',[0.5 0.5 0.5],LW{:},HV{:})
        plot([0 0],[-100 100],'--','Color',[0.5 0.5 0.5],LW{:},HV{:})
        colororder(newcolours)
        set(ax{K,K2},'ColorOrderIndex',1)
        
        switch K2
            case 1
                plot(vktheta0,vs_shear,'Color',map_BrBG(end-1,:),LW{:},'DisplayName','$s_\mathrm{shear}$')
                plot(vktheta0,vs_react,'Color',map_BrBG(2,:),LW{:},'DisplayName','$s_\mathrm{reaction}$')
                plot(vktheta0,vs_react+vs_shear,'Color','k',LW{:},'DisplayName','$s_\mathrm{total}$')
                plot(vktheta0(I),vs_react(I)+vs_shear(I),'x k',MS{:},LW{:},HV{:})
            case 2
                set(ax{K,K2},'ColorOrderIndex',3)
                plot(vktheta0,vktheta,'Color','k',LW{:})
                plot(vktheta0(I),vktheta(I),'x k',MS{:},LW{:},HV{:})
        end
        A.FaceColor=newcolours{1};
        A.LineStyle='none';
        A.FaceAlpha=0.5;
        A.HandleVisibility='off'
    end
    
    vkx0(I)
    vkz0(I)
end
t1.TileSpacing = 'compact';
for K=1:numel(ax)
    ax{K}.TickLabelInterpreter=TX{2};
    set(ax{K},'XLim',[-pi/2 pi/2])
    set(ax{K},FSms{:})
    set(ax{K},LW{:})
    set(ax{K},'XTick',[-pi/2, -3*pi/8 -pi/4, -pi/8, 0, pi/8 pi/4, 3*pi/8, pi/2])
    if K<=2
        set(ax{K},'YTick',[- 20 :10:40],'YLim',[- 20 40]);
        ylab=ylabel(ax{K},'$s$',FSms{:},TX{:},'Rotation',0);
    else
        set(ax{K},'YTick',[-3*pi/4, -pi/2, -pi/4, 0,pi/4],'YLim',[- 3*pi/4 1.3*pi/4]);
        set(ax{K},'YTickLabel',{'$-3\pi/4$','$-\pi/2$','$-\pi/4$','$0$','$\pi/4$'});
        ylab=ylabel(ax{K},'$\theta_1$',FSms{:},TX{:},'Rotation',0);
    end
    xlab=xlabel(ax{K},'Initial angle: $\theta_0$',FSms{:},TX{:});
    set(ax{K},'Layer','top')
    %  set(ax{K},'XTickLabel',{'$-\pi/2$', '$-3\pi/8$','$-\pi/4$','$-\pi/8$','$0$','$\pi/8$','$\pi/4$', '$3\pi/8$', '$\pi/2$'});
    set(ax{K},'XTickLabel',{'$-\pi/2$', [],'$-\pi/4$',[],'$0$',[],'$\pi/4$', [], '$\pi/2$'});
end


l=legend(ax{1});
l.FontSize=FSms{2};
l.Interpreter=TX{2};
l.ItemTokenSize(1)=10;
text(ax{1},-1.5,35,'(a): Near ridge axis',FSms{:},TX{:});
text(ax{2},-1.5,35,'(b): Off ridge axis',FSms{:},TX{:});

f1.Units='centimeters';
f1.Position=[0 0 20 10];
t1.TileSpacing='compact';
t1.Padding='compact';
print(strcat(path_spec,'MOR_growth_spectrum'),output_format)