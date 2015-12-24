function [dispersionAll dispersionClus distances allResults allResultsSD subjMiss rmseMod]=behavAnal2()

%% Behavioral analyses
% Running model analyses on within environment.
close all
clear all
set(0, 'DefaultAxesFontSize',20)
set(0, 'DefaultLineLineWidth',2)
% For mh fits
numBurn2=1000;
numSamps2=200;

labels={'4C1','2C2','1C4','8C1','4C2','2C4','1C8'};
colors={'r','g','b','c','m','y','k'};
% For when >1 object per cluster
labels2={'2C2','1C4','4C2','2C4','1C8'};
colors2={'g','b','m','y','k'};
% For when >1 cluster
labels3={'4C1','2C2','8C1','4C2','2C4'};
colors3={'r','g','c','m','y'}; 


%% Load data
colors={'r','g','b','c','m','y','k'};
colors2={'r','g','c','m','y'};
xcent=500;
ycent=350;
numClu=[4 2 1 8 4 2 1];
numIns=[1 2 4 1 2 4 8];
numItems=[4 4 4 8 8 8 8];
numEnv=10;

loadVarClus;
targets=targs;
%% Set test params
fitDir='behavAnal4_fold';
% behavAnal4_fold: numBurn2=800, numSamps=2400 replicating
% behavAnal3_fold: numBurn2=800, numSamps=2400
% behavAnal2_fold: numBurn2=200, numSamps=600

if ~exist(fitDir)
    mkdir(fitDir)
end
structure={[1 1 1 1],[1 1 2 2],[1 1 1 1],[1 1 1 1 1 1 1 1],[1 1 1 1 2 2 2 2],[1 1 2 2 3 3 4 4],[1 1 1 1 1 1 1 1]};
numSets=size(guesses,1);

%% Run error model
allResults=nan(9,70);
allResultsSD=nan(9,70);
allAssign=cell(35,70);
subjMiss=nan(35,70);
count=1;
numSets=35; % for debugging
priNums=[2 3 5 7 11 13 17 19];
for i=(1:length(numClu))  
    currC=numClu(i);
    currI=numIns(i);
    currItems=numItems(i);
    
    corrClus=reshape(repmat(1:currC,currI,1),currItems,1)';
    inds=(10*(i-1)+1):(10*i);
    for ei=1:numEnv
        if ~exist(fullfile(fitDir,strcat('subjData','_',num2str((i)),'_',num2str((ei)),'.mat')))
            if ei==1
               disp(strcat('Structure',num2str(i))) 
            end
            disp(ei)
            targs=nan(currItems,2,numSets);
            guess=nan(currItems,2,numSets);
            for n=(1:numSets)
                targs(:,:,n)=allStore{inds(ei),1}{1}';
                guess(:,:,n)=guesses{n,inds(ei)};
            end
            [sdGlobe sdClus sdItem weiClus weiItem assignments accept probTarg llk]=mhHierBayes9_v8(targs,guess,currC,numBurn2,numSamps2);
            [pClus pItem]=extractMiss(assignments,corrClus);
            resData=[sdGlobe(numBurn2+1:end) sdClus(numBurn2+1:end) sdItem(numBurn2+1:end) pClus(numBurn2+1:end) pItem(numBurn2+1:end) weiClus(numBurn2+1:end) weiItem(numBurn2+1:end) probTarg(numBurn2+1:end)];
            save(fullfile(fitDir,strcat('subjData','_',num2str((i)),'_',num2str((ei)),'.mat')),'resData','assignments')
        else
            load(fullfile(fitDir,strcat('subjData','_',num2str((i)),'_',num2str((ei)),'.mat')))
        end
        
        for n=1:numSets
            % Find mode assignment
            tempAss=squeeze(assignments(n,:,:));
            tempGsum=sum(repmat(priNums(1:size(tempAss,1))',1,size(tempAss,2)).^tempAss,1);
            mo=find(tempGsum==mode(tempGsum));
            
            allAssign{n,inds(ei)}=tempAss(:,mo(1));
        end
        
        resData2=resData;
        if i==1 || i==4
            % If singleton clusters, combine item and cluster noise
            resData2(:,2)=sqrt((resData(:,2).^2)+(resData(:,3).^2));
            resData2(:,3)=sqrt((resData(:,2).^2)+(resData(:,3).^2));
        elseif i==3 || i==7
            % If single clusters, combine global and cluster noise
            resData2(:,2)=sqrt((resData(:,2).^2)+(resData(:,1).^2));
            resData2(:,1)=sqrt((resData(:,2).^2)+(resData(:,1).^2));
            
        end
        temp=sum(assignments==repmat(1:size(assignments,2),[size(assignments,1),1,size(assignments,3)]),2)/size(assignments,2);
        subjMiss(:,count)=mean(squeeze(temp),2);
        pMiss=1-mean(temp(:));
        allResults(:,count)=[mean(resData2)';pMiss];
        allResultsSD(:,count)=[std(resData2)';mean(std(temp,0,3))];
        count=count+1;
    end
end
% sdg sdc sdi prc pri wec wei pt

%% RMSE

subjDist=squeeze(mean(reshape(distances,35,10,7),2));

rmseMod=rmseSansMiss(targs,guesses,allAssign);
% figure;set(gcf,'color','white');hold on;
% bar(1:7,mean(subjDist));ylabel('RMSE (px)')
% errorbar(1:7,mean(subjDist),std(subjDist)/sqrt(numSets),'k.');hold off; 
% set(gca,'XTickLabel',labels,'XTick',1:length(labels))
%% Were clusters represented relative to the global center?
% If clusters were represented relative to the global center, we would
% expect the dispersion of the objects to be correlated with cluster noise
dispersion=dispEnv(targets,numClu);
dispersionAll=squeeze(dispersion(1,:,1));
dispersionClus=squeeze(dispersion(1,:,2));

yesClusDisp=dispersionClus(dispersionClus~=0);
yesClusClusNoise=allResults(2,dispersionClus~=0);
yesClusClusNoiseSD=allResultsSD(2,dispersionClus~=0);

disp('Reg Corr')
for ci=1:length(yesClusDisp)/10
    %errorbar((yesClusDisp(((ci-1)*10)+1:ci*10)),yesClusClusNoise(((ci-1)*10)+1:ci*10),yesClusClusNoiseSD(((ci-1)*10)+1:ci*10),strcat(colors2{ci},'.'))
    [a b]=corr((yesClusDisp(((ci-1)*10)+1:ci*10))',yesClusClusNoise(((ci-1)*10)+1:ci*10)');
    disp(strcat('r: ', num2str(a),', p: ',num2str(b)))
end
% xlabel('Cluster SD (px)')
% ylabel('Cluster Noise (px)')
% legend('4C1','2C2','8C1','4C2','2C4')
% hold off;

dispersion2=dispEnv2(targets,numClu);
dispersionAll2=squeeze(dispersion2(1,:,1));
dispersionClus2=squeeze(dispersion2(1,:,2));




