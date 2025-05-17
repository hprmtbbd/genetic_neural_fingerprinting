function h = yyplot_matrix(dat,xtL,ytL,yytL,xlab,ylab,yylab,tit_lab,clims,cmap)

if any(isnan(dat(:)))
    cmap(1,:) = [0,0,0];
end

nY = length(ytL); nX = length(xtL);
yt = 1/2:nY-1/2; xt = 1/2:nX-1/2;

colororder({'k','k'})

yyaxis left
imagesc(1/2,1/2,dat,clims)
colormap(cmap)
yticks(yt),yticklabels(ytL)
xticks(xt),set(gca,'xticklabel',xtL,'xticklabelrotation',30)
set(gca,'FontSize',10)
xlabel(xlab,'FontSize',14)
ylabel(ylab,'FontSize',14)

hold on
for ip = 0:2:nY
    plot(0:nX,ip*ones(nX+1,1),'k-','LineWidth',2)
end

for ip = 0:1:nX
    plot(ip*ones(nY+1,1),0:nY,'k-','LineWidth',2)
end

for iy = 1:size(dat,1), for ix = 1:size(dat,2)
    if ~isnan(dat(iy,ix))
        tx = text(xt(ix)-1/6,yt(iy),sprintf('%0.2f',dat(iy,ix)));
        if abs(dat(iy,ix)) > clims(2)*0.6, set(tx,'Color','w'), end
    else
        tx = text(xt(ix)-1/6,yt(iy),'N.S.');
        set(tx,'Color','w')
    end
end, end

hold off

yyaxis right
yy = gca;
set(yy,'YDir','reverse')
ylim([0,nY]),yticks(1:2:nY),yticklabels(strcat({'\fontsize{14} '},yytL)) 
ylabel(yylab,'FontSize',14)

tit = title(tit_lab); tit_pos = get(tit,'Position');
set(tit,'FontSize',14), set(tit,'Position',[tit_pos(1) tit_pos(2)-0.15 tit_pos(3)])

end
