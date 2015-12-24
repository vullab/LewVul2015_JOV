function [fAssign allSd allClusPos]=run_clustering3(currTarg,alpha,sdprior)

%% run_clustering
% Get dirichlet clustering
% 12/14/2014-Saving for different clustering structures
% 12/14/2014-Saving for individual environments

disp('Generating clusters')
numSubj=1;
allZ=cell(numSubj,size(currTarg,2));
allSd=cell(numSubj,size(currTarg,2));
allSdZinds=cell(numSubj,size(currTarg,2));
allEstPos=cell(numSubj,size(currTarg,2));
allClusPos=cell(numSubj,size(currTarg,2));
allLlk=cell(numSubj,size(currTarg,2));
fAssign=cell(numSubj,size(currTarg,2));



for ei=1:70
    dirichName=strcat('dirichSave',num2str(ei),'HG.mat');
    if ~exist(dirichName)           
        disp(ei)
        
        llk=-inf;
        zStore=nan;
        sdStore=nan;
        cluStore=nan;
        
        numChain=40;
        for ic=1:numChain
            [zStore2 sdStore2 llk2 cluStore2]=dirichGibbs9_v3(currTarg{ei},300,100,alpha,sdprior);
            if llk2>llk
                llk=llk2;
                zStore=zStore2;
                sdStore=sdStore2;
                cluStore=cluStore2;
            end
        end
        
        plotClus(currTarg{ei},zStore,cluStore,sdStore);
        title(llk)
        fname=strcat('Env_',num2str(ei));
        saveas(gcf, fname, 'pdf')
        
        
        close
        save(dirichName,'sdStore','zStore','llk','cluStore')
    else
        load(dirichName)
    end
    
    allSd{ei}=sdStore; % sd of the clusters
    allClusPos{ei}=cluStore;
    fAssign{ei}=zStore;
    allLlk{ei}=llk;
end















