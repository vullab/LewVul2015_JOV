%% runDispClus
% Find correlation between cluster dispersion and relative distance

%% Calculate distance between clusters

clusDisp=nan(length(subjs),length(allptime),10);
clusNoise=nan(length(subjs),length(allptime),10);
for is=1:length(subjs)
    for ip=1:length(allptime)
        for ie=1:10
            curr_xt=squeeze(xt(is,ip,ie,:));
            curr_yt=squeeze(yt(is,ip,ie,:));
            
            % Find cluster centers
            clusDisp(is,ip,ie)=sqrt(((mean(curr_xt(1:4))-mean(curr_xt(5:8))).^2)+((mean(curr_yt(1:4))-mean(curr_yt(5:8))).^2));
            
            
            
            clusNoise(is,ip,ie)=allMod2{is,2,ip,ie}(numBurn:end);
        end
    end
end