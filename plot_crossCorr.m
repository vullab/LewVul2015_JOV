%% Run cross correlation 2
% Using product divided by product of norms

for e=1:7
    numO=numC(e)*numI(e);
    temp=zeros(numO,numO);
    for si=1:size(subjCC2,1)
        for ei=1:10
            temp=temp+subjCC2{si,((e-1)*10)+ei};
        end
    end
    temp=temp/(size(subjCC2,1)*10);
    temp=flipud(temp);
    figure('Position', [100, 100, 800, 800]);set(gcf,'Color','white');
    hold on
    imagesc(temp,[-1 1])
    for xi=1:numO
        for yi=1:numO
            if xi+yi==numO+1
                rectangle('Position',[xi-.5,yi-.5,1,1],'FaceColor','k')
            end
        end
    end
    xlim([.5 numO+.5])
    ylim([.5 numO+.5])
    if e~=1 && e~=4
        if e~=3 && e~=7
            for ci=1:(numC(e))
                rectangle('Position',[((ci)-1)*numI(e)+.5,(-numI(e)+1)+numO-((ci)-1)*numI(e)-.5,numI(e),numI(e)],'LineWidth',12,'LineStyle','-','EdgeColor',[0.5  0    0.9])
            end
        else
            for ci=1:(numC(e))
                rectangle('Position',[((ci)-1)*numI(e)+.5+.04,(-numI(e)+1)+numO-((ci)-1)*numI(e)-.5+.04,numI(e)-.09,numI(e)-.09],'LineWidth',12,'LineStyle','-','EdgeColor',[0.5  0    0.9])
            end
        end
    end
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    hold off
end

disp('Were objects more correlated than expected by chance?')
for ei=1:7
    [h p ci stat]=ttest(subjClus(:,ei));
    disp(strcat('t: ',num2str(stat.tstat), ' p: ',num2str(p)) )
end
disp('How about non-clustered objects?')
for ei=1:7
    [h p ci stat]=ttest(subjClus2(:,ei));
    disp(strcat('t: ',num2str(stat.tstat), ' p: ',num2str(p)) )
end
disp('Were clustered objects more correlated?')
for ei=1:7
    [h p ci stat]=ttest(subjClus(:,ei)-subjClus2(:,ei));
    disp(strcat('t: ',num2str(stat.tstat), ' p: ',num2str(p)) )
end
