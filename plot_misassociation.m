%% Object to cluster bias

indCheck=[1 2 3 4 5 6 7];
subjMisassociation=reshape(sum(subjResults(4:5,:)),10,7);

figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on;
errorbar(1:3,mean(subjMisassociation(:,1:3)),std(subjMisassociation(:,1:3))/sqrt(size(subjMisassociation,1)),'k','LineWidth',6);
errorbar(4:7,mean(subjMisassociation(:,4:7)),std(subjMisassociation(:,4:7))/sqrt(size(subjMisassociation,1)),'k','LineWidth',6);
plot([3.5 3.5],[0 1],'r','LineWidth',6)
set(gca,'XTickLabel',labels(indCheck),'XTick',1:length(labels(indCheck)))
hYLabel=ylabel('Misassociation Rate');
ylim([0 .2])

set( gca                       , ...
    'FontName'   , 'Helvetica','FontSize',40 );

set([ hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([ hYLabel]  , ...
    'FontSize'   , 50          );
set([ hYLabel], ...
    'FontName'   , 'Helvetica','Color',[.3 .3 .3] );
set([ hYLabel]  , ...
    'FontSize'   , 60          );
set(gca,'box','off')
set([ hYLabel]  , ...
    'FontSize'   , 50          );
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
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

disp('Did object bias vary')
loadInd=[(ones(size(subjMisassociation,1),3)) 2*(ones(size(subjMisassociation,1),4))];
envInd=repmat(1:7,size(subjMisassociation,1),1);
loadInd2=loadInd(:);
envInd2=envInd(:);

tempTable1=table(subjMisassociation(:),loadInd2,envInd2,'VariableNames',{'miss','load','env'});
lme1 = fitlme(tempTable1,'miss~load*env');
disp(lme1)

% nest=[0 0 ; 1  0 ];
% [p,t,stats,terms]=anovan(subjOCbias(:),{loadInd2 envInd2 },'nested',nest,'display','off');
% disp(strcat('F(',num2str(t{2,3}),',',num2str(t{3,3}),'): ',num2str(t{3,6}),' p:',num2str(t{3,7})))
