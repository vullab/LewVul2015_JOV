function plotClus(locs,zStore,cluStore,sdStore)

numClu=length(unique(zStore));
pCol=linspecer(numClu);

figure('Position', [100, 100, 500, 350]);set(gcf,'color','w');
hold on;
for ic=1:numClu
    sel=zStore==ic;
    plot(locs(sel,1),locs(sel,2),'k.','Color',pCol(ic,:))
    circle(cluStore(1,ic),cluStore(2,ic),sdStore(ic),pCol(ic,:))
    
end
xlim([0 1000])
ylim([0 700])
hold off

    function circle(x,y,r,col)
        %x and y are the coordinates of the center of the circle
        %r is the radius of the circle
        %0.01 is the angle step, bigger values will draw the circle faster but
        %you might notice imperfections (not very smooth)
        ang=0:0.01:2*pi;
        xp=r*cos(ang);
        yp=r*sin(ang);
        plot(x+xp,y+yp,'k-','Color',col);
    end

end