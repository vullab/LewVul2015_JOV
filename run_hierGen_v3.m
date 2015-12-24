function [samps rmse]=run_hierGen_v3(firstTarg,clusAssign,numSamps,objNoise,fname)
%% modFit2_hierGen
% Model fitting of the hierarchical generative model meant to complement behavAnal2
if ~exist(fname)
    numEnv=length(objNoise);
    
    samps=cell(numEnv,numSamps);
    rmse=nan(numEnv,numSamps);
    for ei=1:numEnv
        ei2=ceil(ei/10);
        inds=(1+(10*(ei2-1)):ei2*10);
        encNoise=mean(objNoise(inds));

        for si=1:numSamps
            
            samps{ei,si}=takeSamp(firstTarg{ei},clusAssign{ei},encNoise);
            rmse(ei,si)=mean(sqrt(sum(((samps{ei,si}-firstTarg{ei}).^2),2)));
        end
    end
    save(fname,'samps','rmse')
else
    load(fname);
end

    function estPos=takeSamp(objLocs,clusAssign,objNoise)
        estPos=nan(size(objLocs));
        
        % Find cluster and global noise
        numC=length(unique(clusAssign));
        estClus=nan(numC,2);
        
        clusNoise=nan(1,numC);
        numIn=nan(1,numC);
        clusLoc=nan(2,numC);
        for ci=1:numC
            numIn(ci)=sum(clusAssign==ci);
            clusNoise(ci)=objNoise/sqrt(numIn(ci));
            clusLoc(:,ci)=mean(objLocs(clusAssign==ci,:),1)';
        end
        
        gmean=mean(objLocs,1)';
        gloNoise=sqrt(sum((numIn.*clusNoise.^2))/sum(numIn));
        
        for ci=1:numC
            % Estimate cluster locations
            meanClusPos=((gmean./(gloNoise^2))+(numIn(ci)*clusLoc(:,ci)./ ...
                (clusNoise(ci)^2)))./((1/(gloNoise^2))+(numIn(ci)/(clusNoise(ci)^2)));
            sdClusPos=sqrt(1/(1./(gloNoise^2)+    (numIn(ci)./(clusNoise(ci)^2))    ));
            
            estClus(ci,:)=normrnd(meanClusPos,sdClusPos);
        end
        
        % Estimate object locations
        for ii=1:size(objLocs,1)
            currSel=clusAssign(ii);
            if sum(currSel==clusAssign)==1
                iCurr=objLocs(ii,:)';
                
                % Estimate object locations
                meanEstPos=((iCurr./(objNoise^2))+(estClus(currSel,:)'./ ...
                    (clusNoise(currSel)^2)))./((1/(objNoise^2))+(1/(clusNoise(currSel)^2)));
                sdEstPos=sqrt(1/(1./(objNoise^2)+    (1./(clusNoise(currSel)^2))    ));
                
                estPos(ii,:)=normrnd(meanEstPos,sdEstPos);
            else
                estPos(ii,:)=estClus(currSel,:);
            end
        end
        % figure;hold on;plot(objLocs(:,1),objLocs(:,2),'k.');plot(estPos(:,1),estPos(:,2),'r.');hold off
    end
end

