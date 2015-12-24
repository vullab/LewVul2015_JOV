function [zStore sdStore llk estPos cluStore ]=dirichGibbs9(itemGuess,numSamps,numBurn,alpha,encNoise)
% v7-Generating object samples too, adding object encoding noise
% v8-Introducing cluster bias towards global center
% v9-Turns out I should be using precision, not dispersion

numItems=size(itemGuess,1);
totalSamps=numBurn+numSamps;
llk=nan(totalSamps,1);

gMean=mean(itemGuess)';

% Cluster sd initialization
sdCurr=60;
sdStore=nan(totalSamps,1);
sdStore(1)=sdCurr;

% Cluster assignments
zStore=nan(totalSamps,numItems);
zStore(1,:)=1;
clusParams=(min(itemGuess)+rand()*( max(itemGuess)- min(itemGuess)))';
cluStore=cell(totalSamps,1);
cluStore{1}=clusParams;

% Estimated object positions
estPos=nan(size(itemGuess,1),size(itemGuess,2),totalSamps);
estPos(:,:,1)=(repmat(mean(itemGuess),size(itemGuess,1),1)+itemGuess)/2; % Start off half way between objects and grand mean
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
    
    
    
    % Calculate cluster mean precision
    clusPrec=nan(1,K);
    for j=1:K
       clusPrec(j)=encNoise/sqrt(sum(zStore(i,:)==j)); % Pool object variance
    end
    
    % Calculate global mean precision
    gPrec=sqrt(sum(clusPrec.^2)/length(clusPrec)); % Pool cluster variance
    
    % Update means of clusters
    distsI=[];
    for j=1:K
        % Try changing below line to calculating clusParams to actual
        tempCC=mean(itemGuess(zStore(i,:)==j,:),1)';
        numIn=sum(zStore(i,:)==j);
        meanClusPos=((gMean./(gPrec^2))+(numIn*tempCC./(clusPrec(j)^2)))./((1/(gPrec^2))+(numIn/(clusPrec(j)^2)));
        sdClusPos=sqrt(1/(1./(gPrec^2)+    (numIn./(clusPrec(j)^2))    ));
        
        clusParams(:,j)=normrnd(meanClusPos,sdClusPos);
        
        distsI=[distsI estPos(zStore(i,:)==j,:,i-1)'-repmat(clusParams(:,j),1,sum(zStore(i,:)==j))];
    end
    
    cluStore{i}=clusParams;
    
    % Update sdClus
    dists2=abs(reshape(distsI,numel(distsI),1));
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
    
    %% Update object locations
    meanPos=((itemGuess./(encNoise^2))+(clusParams(:,zStore(i,:))'./(repmat(clusPrec(zStore(i,:))',1,2).^2)))./((1/(encNoise^2))+(1./((repmat(clusPrec(zStore(i,:))',1,2)).^2)));
    sdPos=sqrt(1./((1./(encNoise^2))+(1./(repmat(clusPrec(zStore(i,:))',1,2).^2))   ));
    
    %% Sample objects but restrict from overlapping
    tempEstPos=nan(numItems,2);
    sampOrder=randperm(numItems);
    for j=1:numItems
        temp=normrnd(meanPos(sampOrder(j),:),sdPos(sampOrder(j),:));
        
        distcheck=sum(sqrt(sum((tempEstPos-repmat(temp,numItems,1)).^2,2))<-inf);
        sdstart=sdPos(sampOrder(j),:);
        count=1;
        while distcheck>0
            temp=normrnd(meanPos(sampOrder(j),:),sdstart);
            distcheck=sum(sqrt(sum((tempEstPos-repmat(temp,numItems,1)).^2,2))<60);
            sdstart=sdstart*1.25;
            count=count+1;
            if count>300
               %disp('help') 
            end
        end
        tempEstPos(sampOrder(j),:)=temp;
    end

    estPos(:,:,i)=tempEstPos; 
    
end
llk(1:numBurn)=[];
sdStore(1:numBurn)=[];
zStore(1:numBurn,:)=[];
estPos(:,:,1:numBurn)=[];


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