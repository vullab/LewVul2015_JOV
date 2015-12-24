function [samps, rmse]=run_chunk_v2(firstTarg,clusAssign,numSamps,clusMean,clusSd,fname)

if ~exist(fname)
    numEnv=length(clusMean);
    
    samps=cell(numEnv,numSamps);
    rmse=nan(numEnv,numSamps);
    for ei=1:numEnv
        
        
        for si=1:numSamps
            
            samps{ei,si}=takeSamp(clusAssign{ei},clusMean{ei},clusSd{ei});
            rmse(ei,si)=mean(sqrt(sum(((samps{ei,si}'-firstTarg{ei}).^2),2)));
        end
    end
    save(fname,'samps','rmse')
else
    load(fname);
end
    function samp=takeSamp(clusAssign,clusMean,clusSd)
        
        samp=normrnd(clusMean(:,clusAssign),repmat(clusSd(clusAssign),2,1));
        
        
    end



end







