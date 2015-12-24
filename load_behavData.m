%% load_behavData

loadVarClus;
[dispersionAll dispersionClus distances subjResults subjResultsSD subjMiss rmseMod]=behavAnal2();
load('assignments.mat')
subjDists=squeeze(mean(reshape(distances,35,10,7),2));

% Distribution of locations
pos=[];
for ei=1:70
    pos=[pos;targs{1,ei}];
end

%% all pairwise distances
run_pairDist;

%% Cross cluster analysis
run_crossCorr;

%% Landmark analysis
run_landmarkAnalysis;

%% Put pairwise distance and cross cluster analyses into compatible data formats
run_pairDist_CrossCorr;





