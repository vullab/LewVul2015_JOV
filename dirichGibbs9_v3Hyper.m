function [zStore sdStore llk clusParams alphaStore]=dirichGibbs9_v3Hyper(itemGuess,numSamps,numBurn,sdPrior)
% v7-Generating object samples too, adding object encoding noise
% v8-Introducing cluster bias towards global center
% v9-Turns out I should be using precision, not dispersion
% v9_v2-12/6/2014-Making edits, turns out I hsouldn't have been using
% estPos...
% v9_v3-12/6/2014-Removing clusters at the beginning of each iteration,
% adding in multiple sds
% v9_v3Hyper-4/7/2015-Adding hyperparameter for alpha http://web.cse.ohio-state.edu/~kulis/teaching/788_sp12/DP.learnalpha.pdf

numItems=size(itemGuess,1);
totalSamps=numBurn+numSamps;
llk=nan(totalSamps,1);

gMean=mean(itemGuess)';

% Cluster sd initialization
sdStore=60;
sdPrior=sdPrior/3; %sqrt(mean(var(itemGuess)))/3;

% alpha initialization
alphaStore=.1;

% Cluster assignments
zStore=nan(1,numItems);
zStore(1,:)=1;
clusParams=(min(itemGuess)+rand()*( max(itemGuess)- min(itemGuess)))';

% Final
zFinal=nan;
sdFinal=nan;
llkFinal=-inf; 
cluFinal=nan;

% Estimated object positions
K=1;
for i=2:totalSamps
    if i==numBurn
       disp('') 
    end
    
    %% Fit cluster properties
    tau=randperm(numItems);
    alpha=alphaStore(i-1);
    alphaSum=(numItems)+alpha;
    
    for j=1:length(tau)
        tauCurr=tau(j);
        prevClus=zStore(tauCurr);
        
        % Remove cluster assignment of current item
        zStore(tauCurr)=nan;
        if sum(prevClus==zStore)==0
            clusParams(:,prevClus)=[];
            sdStore(prevClus)=[];
            zStore(zStore>prevClus)=zStore(zStore>prevClus)-1;
            K=K-1;
        end
        
        %% Calculate predictive likelihood of each cluster & new cluster
        newClusStar=mean(itemGuess,1)';
        newClustStarSD=sqrt(mean(var(itemGuess)));

        likelihood=prod(normpdf(repmat(itemGuess(tauCurr,:)',1,size(clusParams,2)+1),[clusParams newClusStar] ...
            ,repmat([sdStore newClustStarSD],2,1)),1);

        prior=[sum(repmat(zStore,size(clusParams,2),1)==repmat((1:size(clusParams,2))',1,numItems),2)' alpha]/alphaSum;
        
        % Sample
        posterior=(likelihood+10^-300).*(prior+10^-300);
        try
        zCurr=randsample(1:length(posterior),1,true,posterior);
        catch
           disp('') 
        end
        newClus=false;
        if zCurr==length(posterior)
            newClus=true;
        end
        zStore(tauCurr)=zCurr;

        % Add new cluster if necessary
        if newClus
            clusParams=[clusParams newClusStar];
            sdStore=[sdStore sdPrior];
            K=K+1;
        end
        
    end 

    %% Update means of clusters
    
    for j=1:K
        % Try changing below line to calculating clusParams to actual
        tempCC=mean(itemGuess(zStore==j,:),1)';
        numIn=sum(zStore==j);
        clusParams(:,j)=tempCC;
        
        distsI=[];
        distsI=[distsI itemGuess(zStore==j,:)'-repmat(clusParams(:,j),1,sum(zStore==j))];
        
        % Update sdClus
        dists2=abs(reshape(distsI,numel(distsI),1));
        dists2=[dists2; sdPrior ;sdPrior;sdPrior];
        dfVar = nansum(dists2.^2);
        df = sum(~isnan(dists2));
        p1 = df/2;
        p2 = dfVar/2;
        sdStore(j)=sqrt(1./gamrnd(p1,1./p2));
    end
    
    %% Sample alpha
    
    % calculate x
%     x_temp=betarnd(alpha+1,numItems);
%     pi_weights=(K-1)/(numItems*(-log(x_temp)));
%     pi_1=pi_weights/(1+pi_weights);
%     pi_2=1-pi_1;
%     
%     alpha_new=(pi_1*gamrnd(K,-log(x_temp))+(pi_2)*gamrnd(K-1,-log(x_temp)));
    alpha_new=1/gamrnd(1+K-1,1+exp(1)+log(numItems));
    alphaStore(i)=alpha_new;
    
    %% Calculate llk
    ll=normpdf(itemGuess',clusParams(:,zStore),repmat(sdStore(zStore),2,1));
    pr=[sum(repmat(zStore,size(clusParams,2),1)==repmat((1:size(clusParams,2))',1,numItems),2)' alpha]/alphaSum;
    llk=prod(prod([ll;pr(zStore)]));
    
    if llkFinal<llk
        llkFinal=llk;
        zFinal=zStore;
        sdFinal=sdStore;
        cluFinal=clusParams;
        
    end
end

end