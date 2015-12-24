function [zStore sdStore llk cluStore ]=dirichGibbs9_v2(itemGuess,numSamps,numBurn,alpha)
% v7-Generating object samples too, adding object encoding noise
% v8-Introducing cluster bias towards global center
% v9-Turns out I should be using precision, not dispersion
% v9_v2-12/6/2014-Making edits, turns out I hsouldn't have been using
% estPos...

numItems=size(itemGuess,1);
totalSamps=numBurn+numSamps;
llk=nan(totalSamps,1);

gMean=mean(itemGuess)';

% Cluster sd initialization
sdCurr=60;
sdStore=nan(totalSamps,1);
sdStore(1)=sdCurr;
sdPrior=sqrt(mean(var(itemGuess)))/3;

% Cluster assignments
zStore=nan(totalSamps,numItems);
zStore(1,:)=1;
clusParams=(min(itemGuess)+rand()*( max(itemGuess)- min(itemGuess)))';
cluStore=cell(totalSamps,1);
cluStore{1}=clusParams;

% Estimated object positions
K=1;
for i=2:totalSamps
    if i==numBurn
       disp('') 
    end
    
    %% Fit cluster properties
    tau=randperm(numItems);
    zStore(i,:)=zStore(i-1,:);
    alphaSum=(numItems)+alpha;
    for j=1:length(tau)
        tauCurr=tau(j);
        prevClus=zStore(i,tauCurr);
        
        % Remove cluster assignment of current item
        zStore(i,tauCurr)=nan;
        
        %% Calculate predictive likelihood of each cluster & new cluster
        newClusStar=[rand()*1000;rand()*700];
        %newClusStar=(min(itemGuess)+rand()*( max(itemGuess)- min(itemGuess)))'; % Assume clusters drawn uniformly over range of object locations
        likelihood=prod(normpdf(repmat(itemGuess(tauCurr,:)',1,size(clusParams,2)+1),[clusParams newClusStar],sdCurr),1);
        prior=[sum(repmat(zStore(i,:),size(clusParams,2),1)==repmat((1:size(clusParams,2))',1,numItems),2)' alpha]/alphaSum;
        
        % Sample
        posterior=(likelihood+10^-300).*(prior+10^-300);
        zCurr=randsample(1:length(posterior),1,true,posterior);

        newClus=false;
        if zCurr==length(posterior)
            newClus=true;
        end
        zStore(i,tauCurr)=zCurr;
        
        % If cluster empty, remove
        if sum(zStore(i,:)==prevClus)==0 && prevClus~=zCurr
            clusParams(:,prevClus)=[];
            zStore(i,:)=shiftDown(zStore(i,:));
            K=K-1;
        end
        
        % Add new cluster if necessary
        if newClus
            clusParams=[clusParams newClusStar];
            K=K+1;
        end
        
    end 

    % Update means of clusters
    distsI=[];
    for j=1:K
        % Try changing below line to calculating clusParams to actual
        tempCC=mean(itemGuess(zStore(i,:)==j,:),1)';
        numIn=sum(zStore(i,:)==j);
        clusParams(:,j)=tempCC;  
        distsI=[distsI itemGuess(zStore(i,:)==j,:)'-repmat(clusParams(:,j),1,sum(zStore(i,:)==j))];
    end
    
    cluStore{i}=clusParams;
    
    % Update sdClus
    dists2=abs(reshape(distsI,numel(distsI),1));
    dists2=[dists2; sdPrior];
    dfVar = nansum(dists2.^2);
    df = sum(~isnan(dists2));
    p1 = df/2;
    p2 = dfVar/2;
    sdCurr=sqrt(1./gamrnd(p1,1./p2));
    sdStore(i)=sdCurr;
    
    % Calculate llk
    ll=normpdf(itemGuess',clusParams(:,zStore(i,:)),sdStore(i));
    pr=sum(repmat(zStore(i,:),size(clusParams,2),1)==repmat((1:size(clusParams,2))',1,numItems),2)'/alphaSum;
    llk(i)=prod([ll(:);pr(:)]);
end
llk(1:numBurn)=[];
sdStore(1:numBurn)=[];
cluStore(1:numBurn)=[];
zStore(1:numBurn,:)=[];



    function raw=shiftDown(raw)
        % Shift z indices so compact
        currInds=unique(raw);
        currInds(isnan(currInds))=[];
        for a=1:length(currInds)
            if raw~=a
                raw(min(raw(raw>a))==raw)=a;
            end
        end
    end

end