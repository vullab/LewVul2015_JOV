%% Find pairwise distance between all objects

if ~exist(fullfile(strcat('subjPairDist','.mat')))
    subjPairDist=cell(length(subjid),size(distances,2));
    for i=1:length(subjid)
        for j=1:size(distances,2)
            currGuess=guesses{i,j};
            currTarg=targs{i,j};
            pdist=findPairDist(currTarg);
            subjPairDist{i,j}=pdist;
            
        end
    end
    save(fullfile(strcat('subjPairDist','.mat')),'subjPairDist')
else
    load(fullfile(strcat('subjPairDist','.mat')))
end