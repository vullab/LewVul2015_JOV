%% main8
% 7/30/2014
% Cleaning up the results from main5 so publishable
% 6_v2-Testing out hybrid model
% 12/6/2014-main8-Recleaning code, using new models. Use this, not main9.
% 2/19/2015-Cleaning for JOV

close all
clear all
clc

% Load settings
labels={'4C1','2C2','1C4','8C1','4C2','2C4','1C8'};
pCol=linspecer(7); % brg

% Environment setting
xcent=500;
ycent=350;
numClu=[4 2 1 8 4 2 1];
numIns=[1 2 4 1 2 4 8];
numItems=[4 4 4 8 8 8 8];
numEnv=10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run behavioral analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running behavioral analyses')

load_behavData;

plot_behavData;

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run model fits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Running model fits')
firstTarg=targs(1,:); % Basic envs

%% Run clustering algorithm on data
disp('Clustering items')
allvar=[];
for ei=1:size(firstTarg,2)
    allvar=[allvar var(firstTarg{ei},[],1)];
end
sdprior=sqrt(mean(allvar));
%[clusAssign clusSd clusMean]=run_clustering3(firstTarg,.1,sdprior);

[clusAssign clusSd clusMean]=run_clustering4(firstTarg,sdprior);
%% Cluster correlation for inferred clusters
plot_crossCorrMod;

%% Model settings
modFold={'modFitsRelPol','modFitsHierGen','modFitsChunk'};
version='_v10';
% _v3-Trying alpha=1, 200 iterations, 10 chains
% _v4-Trying alpha=.01, 100 iterations
% _v5-Trying alpha=.1, 300 iterations, 40 chains
% Trying to match to rmseMod now
% _v6-Trying alpha=1, 300 iterations, 20 chains-Good, but E5 is weird
% _v7-Trying alpha=.5, 300 iterations, 40 chains
% _v8-Trying alpha=.01, 300 iterations, 40 chains-This seems bad
% _v9-Trying alpha=.1, 300 iterations, 20 chains
% _v10-alpha=.1, 300 iterations, 30 chains

numModels=length(modFold);
for mi=1:numModels
    modFold{mi}=strcat(modFold{mi},version);
    if ~exist(modFold{mi})
        mkdir(modFold{mi});
    end
end

%% Run relative position model
% Fit model
relPosFitname=strcat(modFold{1},'_fits','.mat');
relPosFitFname=fullfile(modFold{1},relPosFitname);
[params]=run_fitrelPos(firstTarg,guesses,subjAssign,clusAssign,relPosFitFname);

relPosSampname=strcat(modFold{1},'_samps','.mat');
relPosSampFname=fullfile(modFold{1},relPosSampname);
numSamps=1000;
[samps relpos_rmse]=run_relpos_v1(firstTarg,clusAssign,numSamps,params,relPosSampFname);


%% Run hierarchical generative model

% Generate samples
hierGenSampname=strcat(modFold{2},'_samps','.mat');
hierGenSampFname=fullfile(modFold{2},hierGenSampname);
objNoise=subjResults(3,:);
numSamps=1000;
[samps hiergen_rmse]=run_hierGen_v2(firstTarg,clusAssign,numSamps,objNoise,hierGenSampFname,clusSd);
%[samps hiergen_rmse]=run_hierGen_v2(firstTarg,clusAssign,numSamps,mean(distances,1),hierGenSampFname,clusSd);

% Find rmse
% figure;
% errorbar(mean(hiergen_rmse,2),mean(distances),std(distances)/sqrt(size(distances,1)),'k.')

%% Run chunking model

% Generate samples
chunkSampname=strcat(modFold{3},'_samps','.mat');
chunkSampFname=fullfile(modFold{3},hierGenSampname);
numSamps=1000;
[samps chunk_rmse]=run_chunk_v2(firstTarg,clusAssign,numSamps,clusMean,clusSd,chunkSampFname);



%% Plot model results
%plot_modelFits(rmseMod,relpos_rmse',hiergen_rmse',pCol)
[withCorrs acrCorrs]=plot_modelFits(distances,relpos_rmse',hiergen_rmse',chunk_rmse',pCol);
% figure;plot(mean(relpos_rmse(41:50,:),2),mean(rmseMod(:,41:50)),'k.')






