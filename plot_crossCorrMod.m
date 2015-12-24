%% Run cross correlation 2
% Using product divided by product of norms

subjCorr=nan(35,7,10);
modCorr=nan(35,7,10);
clusAssign;

infClus=[];
for e=1:7
    numO=numC(e)*numI(e);
    uniInd=reshape(repmat(1:numC(e),numI(e),1),1,numO);
    selMat=repmat(uniInd,numO,1)==repmat(uniInd',1,numO);
    if length(unique(uniInd))==numO
        selMat=logical(ones(numO,numO));
    end
    
    selMat(logical(eye(numO)))=0;
    
    for si=1:size(subjCC2,1)
        for ei=1:10
            currMod=clusAssign{10*(e-1)+ei};
            selMatMod=repmat(currMod,numO,1)==repmat(currMod',1,numO);
            if length(unique(currMod))==numO
                selMatMod=logical(ones(numO,numO));
            end
            selMatMod(logical(eye(numO)))=0;
            infClus=[infClus length(unique(currMod))];
            
            subjCorr(si,e,ei)=nanmean(subjCC2{si,((e-1)*10)+ei}(selMat));
            modCorr(si,e,ei)=nanmean(subjCC2{si,((e-1)*10)+ei}(selMatMod));
        end
    end
end

subjCorr2=mean(subjCorr,3);
modCorr2=mean(modCorr,3);
difCorr2=mean(modCorr-subjCorr,3);
%% Plot it


figure('Position', [100, 100, 1049, 895]);set(gcf,'color','w');
hold on
a=errorbar((1:3)+.1,mean(subjCorr2(:,1:3)),std(subjCorr2(:,1:3))/sqrt(numSets),'k','LineWidth',6);
errorbar((4:7)+.1,mean(subjCorr2(:,4:7)),std(subjCorr2(:,4:7))/sqrt(numSets),'k','LineWidth',6);

b=errorbar((1:3)-.1,mean(modCorr2(:,1:3)),std(modCorr2(:,1:3))/sqrt(numSets),'Color',pCol(2,:),'LineWidth',6);
errorbar((4:7)-.1,mean(modCorr2(:,4:7)),std(modCorr2(:,4:7))/sqrt(numSets),'Color',pCol(2,:),'LineWidth',6);

c=errorbar((1:3),mean(difCorr2(:,1:3)),std(difCorr2(:,1:3))/sqrt(numSets),'Color',pCol(3,:),'LineWidth',6);
errorbar((4:7),mean(difCorr2(:,4:7)),std(difCorr2(:,4:7))/sqrt(numSets),'Color',pCol(3,:),'LineWidth',6);

plot([3.5 3.5],[-.1 1],'r','LineWidth',6)
ylim([-.1 1])
set(gca,'XTickLabel',labels,'FontSize',60,'XTick',1:length(labels));
set(gca, 'XColor', [.3 .3 .3]);
hYLabel=ylabel('Error Similarity');
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

legend([a b c],'Actual','Inferred','Difference')
hold off;
