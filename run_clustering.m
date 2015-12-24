function [fAssign allSd allClusPos]=run_clustering(currTarg,alpha,sdprior)

%% run_clustering
% Get dirichlet clustering

numSubj=1;
dirichName='dirichSave6HG.mat';
if ~exist(dirichName)
    disp('Generating clusters')
    
    allZ=cell(numSubj,size(currTarg,2));
    allSd=cell(numSubj,size(currTarg,2));
    allSdZinds=cell(numSubj,size(currTarg,2));
    allEstPos=cell(numSubj,size(currTarg,2));
    allClusPos=cell(numSubj,size(currTarg,2));
    allLlk=cell(numSubj,size(currTarg,2));
    fAssign=cell(numSubj,size(currTarg,2));
    
    for si=1:numSubj
        disp(strcat('Subj ',num2str(si)))
        for ei=1:size(currTarg,2)
            disp(ei)
            
            llk=-inf;
            zStore=nan;
            sdStore=nan;
            cluStore=nan;
            
            numChain=20;
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
            
            % disp(zStore)
            allSd{si,ei}=sdStore; % sd of the clusters
            allClusPos{si,ei}=cluStore;
            fAssign{si,ei}=zStore;
            allLlk{si,ei}=max(llk);
            close 
        end
    end
    save(dirichName,'allSd','fAssign','allLlk','allClusPos')
else
    load(dirichName)
end











