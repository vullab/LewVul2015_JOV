function [fAssign allSd allClusPos]=run_clustering4(currTarg,sdprior)

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


allAlpha=nan(1,70);
for ei=1:70
    dirichName=strcat('dirichSave_Hyper',num2str(ei),'HG.mat');
    if ~exist(dirichName)           
        disp(ei)
        
        llk=-inf;
        zStore=nan;
        sdStore=nan;
        cluStore=nan;
        
        numChain=40;
        for ic=1:numChain
            [zStore2 sdStore2 llk2 cluStore2 alphaStore2]=dirichGibbs9_v3Hyper(currTarg{ei},600,100,sdprior);
            if llk2>llk
                llk=llk2;
                zStore=zStore2;
                sdStore=sdStore2;
                cluStore=cluStore2;
                alphaStore=alphaStore2;
            end
        end
        
        plotClus(currTarg{ei},zStore,cluStore,sdStore);
        title(llk)
        fname=strcat('EnvHyper_',num2str(ei));
        saveas(gcf, fname, 'pdf')
        
        
        close
        save(dirichName,'sdStore','zStore','llk','cluStore','alphaStore')
    else
        load(dirichName)
    end
    allSd{ei}=sdStore; % sd of the clusters
    allClusPos{ei}=cluStore;
    fAssign{ei}=zStore;
    allLlk{ei}=llk;
    allAlpha(ei)=median(alphaStore(100:end));
end


allAlpha;













