%% Object to cluster bias

indCheck=[2 3 5 6 7];
subjOCbias=1-reshape(subjResults(7,:),10,7);

figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on;
errorbar(1:2,mean(subjOCbias(:,2:3)),std(subjOCbias(:,2:3))/sqrt(size(subjOCbias,1)),'k','LineWidth',6);
errorbar(3:5,mean(subjOCbias(:,5:7)),std(subjOCbias(:,5:7))/sqrt(size(subjOCbias,1)),'k','LineWidth',6);
plot([2.5 2.5],[0 1],'r','LineWidth',6)
set(gca,'XTickLabel',labels(indCheck),'XTick',1:length(labels(indCheck)))
hYLabel=ylabel('Object-to-Cluster Bias');
ylim([0 .5])

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

% disp('Did cluster to global bias vary')
% subjOCbias=1-reshape(subjResults(6,:),10,7);
% loadInd=[(ones(size(subjOCbias,1),3)) 2*(ones(size(subjOCbias,1),4))];
% envInd=repmat(1:7,size(subjOCbias,1),1);
% subjOCbias(:,[3 7])=[];
% loadInd(:,[3 7])=[];
% envInd(:,[3 7])=[];
% loadInd2=loadInd(:);
% envInd2=envInd(:);
% 
% tempTable1=table(subjOCbias(:),loadInd2,envInd2,'VariableNames',{'bias','load','env'});
% lme1 = fitlme(tempTable1,'bias~load+env');
% disp(lme1)


% nest=[0 0 ; 1  0 ];
% [p,t,stats,terms]=anovan(subjOCbias(:),{loadInd2 envInd2 },'nested',nest,'display','off');
% disp(strcat('F(',num2str(t{2,3}),',',num2str(t{3,3}),'): ',num2str(t{3,6}),' p:',num2str(t{3,7})))

%% Pairwise comparisons
subjOCbias2=subjOCbias(:);

subjOCbias2_4obj=subjOCbias2(1:(20));
envInd2_4obj=envInd2(1:(20));
[p,t,stats,terms]=anovan(subjOCbias2_4obj(:),{envInd2_4obj   },'display','off');
disp('4 object load multiple comparison')
c_4obj=multcompare(stats,'CType','hsd');
disp(c_4obj)

subjOCbias2_8obj=subjOCbias2((21):end);
envInd2_8obj=envInd2((21):end);
[p,t,stats,terms]=anovan(subjOCbias2_8obj(:),{envInd2_8obj   },'display','off');
disp('8 object load multiple comparison')
c_8obj=multcompare(stats,'CType','hsd');
disp(c_8obj)



