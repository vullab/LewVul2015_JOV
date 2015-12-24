%% plot_rmse

figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on
errorbar(1:3,mean(subjDists(:,1:3)),std(subjDists(:,1:3))/sqrt(numSets),'k','LineWidth',6);
errorbar(4:7,mean(subjDists(:,4:7)),std(subjDists(:,4:7))/sqrt(numSets),'k','LineWidth',6);
plot([3.5 3.5],[0 150],'r','LineWidth',6)
ylim([0 150])
set(gca,'XTickLabel',labels,'FontSize',60,'XTick',1:length(labels));
set(gca, 'XColor', [.3 .3 .3]);
hYLabel=ylabel('RMSE (px)');
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

disp('Did RMSE vary with load and structure')
subjInd=repmat((1:size(subjDists,1))',1,size(subjDists,2));
loadInd=[(ones(size(subjDists,1),3)) 2*(ones(size(subjDists,1),4))];
envInd=repmat(1:7,size(subjDists,1),1);
subjInd2=subjInd(:);
loadInd2=loadInd(:);
envInd2=envInd(:);

tempTable1=table(subjDists(:),subjInd2,loadInd2,envInd2,'VariableNames',{'dists','subs','load','env'});
% lme1 = fitlme(tempTable1,'dists~load+env+(1|subs)+(1|subs:env)+(1|subs:load)');
% disp(lme1)

lme1 = fitlme(tempTable1,'dists~load*env+(1|subs)+(1|subs:env)+(1|subs:load)+(1|subs:load:env)');
disp(lme1)

nest=[0 0 0; 1 0 0;0 0 0];
[p,t,stats,terms]=anovan(subjDists(:),{loadInd2 envInd2 subjInd2  },'nested',nest,'display','off');
disp(strcat('F(',num2str(t{2,3}),',',num2str(t{3,3}),',',num2str(t{4,3}),'): ',num2str(t{3,6}),' p:',num2str(t{3,7})))

%% Individual comparisons
subjDists2=subjDists(:);

subjDists2_4obj=subjDists2(1:(3*35));
envInd2_4obj=envInd2(1:(3*35));
subjInd2_4obj=subjInd2(1:(3*35));
[p,t,stats,terms]=anovan(subjDists2_4obj(:),{envInd2_4obj subjInd2_4obj  },'display','off');
disp('4 object load multiple comparison')
c_4obj=multcompare(stats,'CType','hsd');
disp(c_4obj)

subjDists2_8obj=subjDists2(((3*35)+1):end);
envInd2_8obj=envInd2(((3*35)+1):end);
subjInd2_8obj=subjInd2(((3*35)+1):end);
[p,t,stats,terms]=anovan(subjDists2_8obj(:),{envInd2_8obj subjInd2_8obj  },'display','off');
disp('8 object load multiple comparison')
c_8obj=multcompare(stats,'CType','hsd');
disp(c_8obj)


