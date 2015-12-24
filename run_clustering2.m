function [fAssign allSd allClusPos]=run_clustering2(currTarg,alpha,sdprior)

%% run_clustering
% Get dirichlet clustering
% 12/14/2014-Saving for different clustering structures


disp('Generating clusters')
numSubj=1;
allZ=cell(numSubj,size(currTarg,2));
allSd=cell(numSubj,size(currTarg,2));
allSdZinds=cell(numSubj,size(currTarg,2));
allEstPos=cell(numSubj,size(currTarg,2));
allClusPos=cell(numSubj,size(currTarg,2));
allLlk=cell(numSubj,size(currTarg,2));
fAssign=cell(numSubj,size(currTarg,2));

for si=1:7
    dirichName=strcat('dirichSave',num2str(si),'HG.mat');
    rang=((si-1)*10+1):((si-1)*10+10);
    if ~exist(dirichName)
        zTemp=cell(1,10);
        sdTemp=cell(1,10);
        cluTemp=cell(1,10);
        llkTemp=cell(1,10);
        
        for ei2=1:10
            ei=(si-1)*10+ei2;
            disp(ei)
            
            llk=-inf;
            zStore=nan;
            sdStore=nan;
            cluStore=nan;
            
            numChain=60;
            numIts=300;
            for ic=1:numChain
                [zStore2 sdStore2 llk2 cluStore2]=dirichGibbs9_v3(currTarg{ei},numIts,100,alpha,sdprior);
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
            
            zTemp{1,ei2}=zStore;
            sdTemp{1,ei2}=sdStore;
            cluTemp{1,ei2}=cluStore;
            llkTemp{1,ei2}=llk;
            
            close
        end
        save(dirichName,'sdTemp','zTemp','llkTemp','cluTemp','numChain','numIts')
    else
        load(dirichName)
    end
    allSd(1,rang)=sdTemp; % sd of the clusters
    allClusPos(1,rang)=cluTemp;
    fAssign(1,rang)=zTemp;
    allLlk(1,rang)=llkTemp;
end












