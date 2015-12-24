%% Plot landmark results
figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on

%errorbar(mean(clusPosX,2)-500,mean(clusErrX,2),std(clusErrX,[],2)./sqrt(size(clusErrX,2)),'k.')

%% In X
[a b]=hist(mean(clusPosX,2),10);
binDiff=b(2)-b(1);
[a2 b2]=histc(mean(clusPosX,2),[0 b+binDiff]);
histline=nan(1,10);
histlinesem=nan(1,10);
histX=[];
histXid=[];
for i=1:10
    histX=[histX ; mean(clusErrX(b2==i,:),2)];
    histXid=[histXid ; repmat(i,length(mean(clusErrX(b2==i,:),2)),1)];
    histline(i)=mean(mean(clusErrX(b2==i,:),2));
    histlinesem(i)=std(mean(clusErrX(b2==i,:),2))./sqrt(sum(b2==i));
end
disp('Errors in X dimension?')
tempTable1=table(histX,histXid,'VariableNames',{'rmse','bin'});
lme1 = fitlme(tempTable1,'rmse~bin');
disp(lme1)

% [p,t,stats,terms]=anovan(histX,{histXid },'display','off');
% disp(strcat('F(',num2str(t{2,3}),'): ',num2str(t{2,6}),' p:',num2str(t{2,7})))

errorbar(b-500,histline,histlinesem,'.-','Color',[.3 .3 .3],'LineWidth',6)

xlim([-500 500])

%% In Y
[a b]=hist(mean(clusPosY,2),10);
binDiff=b(2)-b(1);
[a2 b2]=histc(mean(clusPosY,2),[0 b+binDiff]);
histline=nan(1,10);
histlinesem=nan(1,10);
histY=[];
histYid=[];
for i=1:10
    histY=[histY ; mean(clusErrY(b2==i,:),2)];
    histYid=[histYid ; repmat(i,length(mean(clusErrY(b2==i,:),2)),1)];
    histline(i)=mean(mean(clusErrY(b2==i,:),2));
    histlinesem(i)=std(mean(clusErrY(b2==i,:),2))./sqrt(sum(b2==i));
end

disp('Errors in Y dimension?')
tempTable1=table(histY,histYid,'VariableNames',{'rmse','bin'});
lme1 = fitlme(tempTable1,'rmse~bin');
disp(lme1)

% disp('Errors in Y dimension?')
% [p,t,stats,terms]=anovan(histY,{histYid },'display','off');
% disp(strcat('F(',num2str(t{2,3}),'): ',num2str(t{2,6}),' p:',num2str(t{2,7})))

errorbar(-1*(b-350),histline,histlinesem,'-','Color',[0 0 0],'LineWidth',6)

%% Details
legend('X-Dim','Y-Dim')
plot([0 0],[40 80],'Color',[.75 .75 .75])

hXLabel=xlabel('Position (px)');
hYLabel=ylabel('Absolute Error (px)');
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
hold off
