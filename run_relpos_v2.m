function [samps rmse]=run_relpos_v2(firstTarg,clusAssign,numSamps,params,fname)

if ~exist(fname)
    numEnv=length(firstTarg);
    
    samps=cell(numEnv,numSamps);
    rmse=nan(numEnv,numSamps);
    for ei=1:numEnv
        currParam=params{ceil(ei/10)}{1};
        
        for si=1:numSamps
            
            samps{ei,si}=takeSamp(firstTarg{ei},clusAssign{ei},currParam);
            rmse(ei,si)=mean(sqrt(sum(((samps{ei,si}-firstTarg{ei}).^2),2)));
        end
    end
    save(fname,'samps','rmse')
else
    load(fname);
end

    function estPos=takeSamp(objLocs,clusAssign,params)
        gmean=mean(objLocs,1);
        numC=length(unique(clusAssign));
        cmean=nan(numC,2);
        for ci=1:numC
            cmean(ci,:)=mean(objLocs(clusAssign==ci,:),1);
        end
        crel=cmean-repmat(gmean,numC,1);
        [cAng cRad]=cart2pol(crel(:,1),crel(:,2));
        
        irel=objLocs-cmean(clusAssign,:);
        [iAng iRad]=cart2pol(irel(:,1),irel(:,2));
        
        for ci=1:numC
            cAng2(ci,1)=circ_vmrnd(cAng(ci),1./(params(2)^2),1);
        end
        for ii=1:size(objLocs,1)
            iAng2(ii,1)=circ_vmrnd(iAng(ii),1./(params(2)^2),1);
        end
        cRad2=exp(normrnd(log(cRad),params(1)));
        iRad2=exp(normrnd(log(iRad),params(1)));
        
        [cx cy]=pol2cart(cAng2,cRad2);
        crel2=[cx cy];
        [ix iy]=pol2cart(iAng2,iRad2);
        irel2=[ix iy];
        
        estPos=repmat(gmean,size(objLocs,1),1)+crel2(clusAssign,:)+irel2;
    end

end