%% Noise of cluster positions given dispersion

indCheck=[1 2     4 5 6];
subjCluNoise=reshape(subjResults(2,:),10,7);

figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
yesClusDisp=dispersionClus;
yesClusClusNoise=subjResults(2,:);
yesClusClusNoiseSD=subjResultsSD(2,:);
hold on;
for ci=1:length(indCheck)
    fill([0 0 0 0],[0 0 0 0],pCol(indCheck(ci),:))
    
end
[legh,objh,outh,outm] =legend('4C1','2C2','8C1','4C2','2C4','Location','NorthWest');

for ci=1:length(indCheck)
    errorbar((yesClusDisp(((indCheck(ci)-1)*10)+1:indCheck(ci)*10)),yesClusClusNoise(((indCheck(ci)-1)*10)+1:indCheck(ci)*10),yesClusClusNoiseSD(((indCheck(ci)-1)*10)+1:indCheck(ci)*10),'k.','Color',pCol(indCheck(ci),:),'LineWidth',6)
    [a b]=corr((yesClusDisp(((indCheck(ci)-1)*10)+1:indCheck(ci)*10))',yesClusClusNoise(((indCheck(ci)-1)*10)+1:indCheck(ci)*10)');
    disp(strcat('r: ',num2str(a),' p: ',num2str(b)))
end
hXLabel=xlabel('Cluster SD (px)');
hYLabel=ylabel('Cluster location noise (px)');


set([legh],'FontName'   , 'Helvetica','TextColor',[.3 .3 .3] )

set([legh],'FontName'   , 'Helvetica','TextColor',[.3 .3 .3] )
set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );
set( [legh]                       , ...
    'FontName'   , 'Helvetica','FontSize',30 );

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

disp('Did cluster dispersion affect cluster accuracy?')
[a b]=corr(yesClusDisp(yesClusDisp~=0)',yesClusClusNoise(yesClusDisp~=0)');
disp(strcat('r= ',num2str(a),' p=',num2str(b)))


disp('Object noise')
subjOCbias=reshape(subjResults(3,:),10,7);
loadInd=[(ones(size(subjOCbias,1),3)) 2*(ones(size(subjOCbias,1),4))];
envInd=repmat(1:7,size(subjOCbias,1),1);
subjOCbias(:,[1 4])=[];
loadInd(:,[1 4])=[];
envInd(:,[1 4])=[];
loadInd2=loadInd(:);
envInd2=envInd(:);

tempTable1=table(subjOCbias(:),loadInd2,envInd2,'VariableNames',{'bias','load','env'});
lme1 = fitlme(tempTable1,'bias~load+env');
disp(lme1)

disp('Cluster noise')
subjOCbias=reshape(subjResults(2,:),10,7);
loadInd=[(ones(size(subjOCbias,1),3)) 2*(ones(size(subjOCbias,1),4))];
envInd=repmat(1:7,size(subjOCbias,1),1);
subjOCbias(:,[3 7])=[];
loadInd(:,[3 7])=[];
envInd(:,[3 7])=[];
loadInd2=loadInd(:);
envInd2=envInd(:);

tempTable1=table(subjOCbias(:),loadInd2,envInd2,'VariableNames',{'bias','load','env'});
lme1 = fitlme(tempTable1,'bias~load+env');
disp(lme1)

disp('Global noise')
subjOCbias=reshape(subjResults(1,:),10,7);
loadInd=[(ones(size(subjOCbias,1),3)) 2*(ones(size(subjOCbias,1),4))];
envInd=repmat(1:7,size(subjOCbias,1),1);
loadInd2=loadInd(:);
envInd2=envInd(:);

tempTable1=table(subjOCbias(:),loadInd2,envInd2,'VariableNames',{'bias','load','env'});
lme1 = fitlme(tempTable1,'bias~load+env');
disp(lme1)


