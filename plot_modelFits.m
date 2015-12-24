function [allCorrs acrCorrs]=plot_modelFits(distances,relpos_rmse,hiergen_rmse,chunk_rmse,pCol)

allCorrs=nan(7,3);
acrCorrs=nan(1,3);

disp('Chunking')
figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on
for ei=1:7
    sel=(ei-1)*10+1:ei*10;
    errorbar(mean(chunk_rmse(:,sel)),mean(distances(:,sel)), ...
        std(distances(:,sel))./sqrt(size(distances(:,sel),1)),'k.','Color',pCol(ei,:),'LineWidth',4);
    [r p]=corr(mean(chunk_rmse(:,sel))',mean(distances(:,sel))');
    disp(strcat('E',num2str(ei),' r: ',num2str(r),'_','p: ',num2str(p)))
    allCorrs(ei,1)=r;
end
plot(0:250,0:250,'k--')
[r p ]=corr(mean(chunk_rmse(:,:))',mean(distances(:,:))');
acrCorrs(1)=r;
[r1 p1 rl1 ru1 ]=corrcoef(mean(chunk_rmse(:,:))',mean(distances(:,:))');
disp(strcat('All',' r: ',num2str(r),'_','p: ',num2str(p)))
disp(strcat('Low',num2str(rl1(1,2)),'_','High',num2str(ru1(1,2))))

hXLabel=xlabel('Chunk RMSE (px)');
hYLabel=ylabel('Subject RMSE (px)');

set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );
set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );

set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 60          );
set(gca,'box','off')
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );

scale = 0.15;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

hold off;

disp('Relative position')
figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on
for ei=1:7
    sel=(ei-1)*10+1:ei*10;
    errorbar(mean(relpos_rmse(:,sel)),mean(distances(:,sel)), ...
        std(distances(:,sel))./sqrt(size(distances(:,sel),1)),'k.','Color',pCol(ei,:),'LineWidth',4);
    [r p]=corr(mean(relpos_rmse(:,sel))',mean(distances(:,sel))');
    disp(strcat('E',num2str(ei),' r: ',num2str(r),'_','p: ',num2str(p)))
    allCorrs(ei,2)=r;
end
plot(0:250,0:250,'k--')
[r p ]=corr(mean(relpos_rmse(:,:))',mean(distances(:,:))');
acrCorrs(2)=r;
[r1 p1 rl1 ru1 ]=corrcoef(mean(relpos_rmse(:,:))',mean(distances(:,:))');
disp(strcat('All',' r: ',num2str(r),'_','p: ',num2str(p)))
disp(strcat('Low',num2str(rl1(1,2)),'_','High',num2str(ru1(1,2))))

hXLabel=xlabel('Rel. Pos. RMSE (px)');
hYLabel=ylabel('Subject RMSE (px)');

set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );
set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );

set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 60          );
set(gca,'box','off')
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );

scale = 0.15;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

hold off;

disp('Hierarchical generative')
figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on

for ei=1:7
    fill([0 0 0 0],[0 0 0 0],pCol(ei,:))
end
legend({'4C1','2C2','1C4','8C1','4C2','2C4','1C8'},'FontSize',36, ...
    'FontName'   , 'Helvetica','TextColor',[.3 .3 .3],'Location','NorthWest')

for ei=1:7
    sel=(ei-1)*10+1:ei*10;
    errorbar(mean(hiergen_rmse(:,sel)),mean(distances(:,sel)), ...
        std(distances(:,sel))./sqrt(size(distances(:,sel),1)),'k.','Color',pCol(ei,:),'LineWidth',4);
    [r p]=corr(mean(hiergen_rmse(:,sel))',mean(distances(:,sel))');
    disp(strcat('E',num2str(ei),' r: ',num2str(r),'_','p: ',num2str(p)))
    allCorrs(ei,3)=r;
end
plot(0:200,0:200,'k--')
xlim([0 200])
ylim([0 200])
[r p]=corr(mean(hiergen_rmse(:,:))',mean(distances(:,:))');
acrCorrs(3)=r;
[r1 p1 rl1 ru1 ]=corrcoef(mean(hiergen_rmse(:,:))',mean(distances(:,:))');
disp(strcat('All',' r: ',num2str(r),'_','p: ',num2str(p)))
disp(strcat('Low',num2str(rl1(1,2)),'_','High',num2str(ru1(1,2))))
hXLabel=xlabel('Hier. Gen. RMSE (px)');
hYLabel=ylabel('Subject RMSE (px)');

set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );
set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );

set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set([hXLabel hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([hXLabel hYLabel]  , ...
    'FontSize'   , 60          );
set(gca,'box','off')
set([hXLabel hYLabel]  , ...
    'FontSize'   , 50          );
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'YGrid'       , 'on'      , ...
    'XColor'      , [.3 .3 .3], ...
    'YColor'      , [.3 .3 .3], ...
    'LineWidth'   , 1         );

scale = 0.15;
pos = get(gca, 'Position');
pos(2) = pos(2)+scale*pos(4);
pos(4) = (1-scale)*pos(4);
set(gca, 'Position', pos)

hold off;


end