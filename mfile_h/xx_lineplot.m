function fig = xx_lineplot(mats,nplot,xt1,xt2,xt3,xt_lab1,xt_lab2,xt_lab3,y_lab,tit_labs,leg_labs,rb_flag)

fig = figure('Position',[20,50,1500,600]);
for iplot = 1:nplot
    mat0 = mats{iplot};
    xlim0 = size(mat0,1);
    
    subplot(nplot/2,2,iplot)
    axUnder = cla();
    axUnder2 = copyobj(axUnder, ancestor(axUnder,'figure'));
    ax = copyobj(axUnder, ancestor(axUnder,'figure'));
    
    h = plot(ax,mat0,'-','LineWidth',1.5);
    % h = plot(ax,mat0,'*-','LineWidth',1.5,'MarkerSize',4);
    if size(mat0,2) == 2
        if ~rb_flag
            set(h(1),'Color',[0.4940 0.1840 0.5560])
            set(h(2),'Color',[0.4660 0.6740 0.1880])
        else
            set(h(1),'Color','r'),set(h(2),'Color','b')
        end
    elseif size(mat0,2) == 3
        set(h(1),'Color',[0.4940 0.1840 0.5560])
        set(h(2),'Color',[0.9290 0.6940 0.1250])
        set(h(3),'Color',[0.4660 0.6740 0.1880])
    end
    set(ax,'XTick',xt1,'XLim',[1,xlim0],'XTickLabel',xt_lab1,'TickDir','in')
    set(ax,'YTick',0:0.1:1,'YLim',[0,1])
    set(ax,'XGrid','on')
    
    %{
    yvec = 0:0.1:1;
    hold(ax,'on')
    for ix = 3:2:(length(xt1)-1)
        plot(ax,xt1(ix)*ones(length(yvec),1),yvec,'k-','LineWidth',1)
    end
    hold(ax,'off')
    %}
    
    ax_pos = get(ax,'Position');
    if ismember(iplot,[1,3]), xpos = ax_pos(1)-0.09;
    elseif ismember(iplot,[2,4]), xpos = ax_pos(1)-0.03;
    end
    if ismember(iplot,[1,2]), ypos = ax_pos(2)+0.07;
    elseif ismember(iplot,[3,4]), ypos = ax_pos(2)+0.12;
    end
    set(ax,'Position',[xpos ypos ax_pos(3)*1.35 ax_pos(4)*0.85])
    
    plot(axUnder, nan, nan, 'x')
    set(axUnder,'YTick',[])
    set(axUnder,'XTick',xt2,'XLim',[1,xlim0],'XTickLabel',xt_lab2)
    set(axUnder,'XTickLabelRotation',0)
    
    plot(axUnder2, nan, nan, 'x')
    set(axUnder2,'YTick',[])
    set(axUnder2,'XTick',xt3,'XLim',[1,xlim0],'XTickLabel',xt_lab3)
    set(axUnder2,'XTickLabelRotation',0)
    axUnder2.XRuler.TickLabelGapOffset = 15;
    
    ylabel(ax,y_lab)
    title(ax,tit_labs{iplot})
    set(ax,'FontSize',12)
    
    set(axUnder,'FontSize',12)
    set(axUnder2,'FontSize',12)
    
    linkprop([ax,axUnder],{'Position','InnerPosition','XLim'});
    linkprop([ax,axUnder2],{'Position','InnerPosition','XLim'});
    axUnder.HandleVisibility = 'off';
    axUnder2.HandleVisibility = 'off';
    
    if iplot == 3
        lgd = legend(leg_labs,'Location','south','orientation','vertical');
        set(lgd,'FontSize',12)
        pos_l = get(lgd,'Position');
        if size(mat0,2) == 2, ypos_l = pos_l(2)*0.2;
        elseif size(mat0,2) == 3, ypos_l = pos_l(2)*0.1;
        end
        new_pos_l = [pos_l(1) ypos_l pos_l(3) pos_l(4)];
        set(lgd,'Position',new_pos_l)
        
        ax3_pos = get(ax,'Position');
    end
end

txt1 = {'Limbic: Insula','             Cingulate'};
txt_box1 = annotation('textbox','String',txt1,'EdgeColor','none','FontSize',12);
set(txt_box1,'FitBoxToText','on')
drawnow
set(txt_box1,'Position',[ax3_pos(1)+0.3045 ax3_pos(2)*0.2+0.085 1 0])

txt2 = {'Subcortical: Mesial Temporal','                    Basal Ganglia',...
    '                    Thalamus'};
txt_box2 = annotation('textbox','String',txt2,'EdgeColor','none','FontSize',12);
set(txt_box2,'FitBoxToText','on')
drawnow
set(txt_box2,'Position',[ax3_pos(1)+0.39 ax3_pos(2)*0.2+0.085 1 0])

end

