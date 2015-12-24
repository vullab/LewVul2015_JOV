function [sdGlobe sdClus sdItem weiClus weiItem assignments accept probTarg llkAll]=mhHierBayes9_v8(targs,guesses,numC,numBurn,numSamps)
% v1-hierarchical bayesian fit
% v2-Adding assignments
% v3-Adjusting for different cluster structures
% v6-Assignment pdf based on num misses
% v7-Adding normalizations term
% v8-Normalization term + prior probs



%% Set up
numI=size(targs,1);
numIn=numI/numC;
numSubj=size(targs,3);
corrClus=reshape(repmat((1:numC)',1,numI/numC)',numI,1);

% Sampler properties
totalSamps=numBurn+numSamps;
llkAll=nan(totalSamps,1);

% Get targ/item rel pos
targCent=repmat(mean(targs),numI,1);
targClus=reshape(repmat(mean(reshape(targs,numIn,2,numC,numSubj),1),numIn,1),numI,2,numSubj)-targCent;
targItem=targs-targClus-targCent;

%% Initialize parameters

% Global sd initialization
sdGlobe=nan(totalSamps,1); sdGlobe(1)=30;

% Cluster sd initialization
sdClus=nan(totalSamps,1); sdClus(1)=30;

% Item sd initialization
sdItem=nan(totalSamps,1); sdItem(1)=30;

% Cluster weight initialization
weiClus=nan(totalSamps,1); weiClus(1)=.75;

% Item weight initialization
weiItem=nan(totalSamps,1); weiItem(1)=.75;

% Prob target initialization
probTarg=nan(totalSamps,1); probTarg(1)=.95;%rand()*.5+.25;

% Assignment initialization
assignments=nan(numSubj,numI,totalSamps);
for i=1:numSubj
    assignments(i,:,1)=randperm(numI);
end
assignments=repmat(1:numI,[numSubj,1,totalSamps]);

%% Run sampler
llkTotal=calcLikelihood4(targCent,targClus,targItem,guesses,sdGlobe(1),sdClus(1),sdItem(1),weiClus(1),weiItem(1),numC,corrClus,assignments(:,:,1),probTarg(1));
llkAll(1)=mean(prod(llkTotal,2));
accept=0;
for i=2:totalSamps
    if i==numBurn
        disp('')
    end
    
    % Get parameters for this iteration
    currSdGlobe=sdGlobe(i-1);
    currSdClus=sdClus(i-1);
    currSdItem=sdItem(i-1);
    currWclus=weiClus(i-1);
    currWitem=weiItem(i-1);
    currProbTarg=probTarg(i-1);
    currAssign=assignments(:,:,i-1);
    
    %% Determine subject params
    
    % Perturbate values
    pShift2=.025;
    dShift=2.5;
    candSdGlobe=bound3(currSdGlobe,dShift);
    candSdClus=bound3(currSdClus,dShift);
    candSdItem=bound3(currSdItem,dShift);
    candWclus=bound(currWclus,.1);%bound(currWclus,pShift2);
    candWitem=bound(currWitem,.1);%bound(currWitem,pShift2);
    candProbTarg=bound(currProbTarg,.1);
    % candSdGlobe=currSdGlobe;candSdClus=currSdClus;candSdItem=currSdItem;candWclus=currWclus;candWitem=0;
    llkCand=calcLikelihood4(targCent,targClus,targItem,guesses,candSdGlobe,candSdClus,candSdItem,candWclus,candWitem,numC,corrClus,currAssign,candProbTarg);

    % Determine if should switch
    a=prod(prod((llkCand+(10^-100))./(llkTotal+(10^-100)),2));
    check=rand()<a;
    if a>1 || check
        if i<numBurn+1
            accept=accept+1;
        end
        sdGlobe(i)=candSdGlobe;
        sdClus(i)=candSdClus;
        sdItem(i)=candSdItem;
        weiClus(i)=candWclus;
        weiItem(i)=candWitem;
        probTarg(i)=candProbTarg;
        llkTotal=llkCand;
        llkAll(i)=mean(prod(llkCand,2));
    else
        sdGlobe(i)=currSdGlobe;
        sdClus(i)=currSdClus;
        sdItem(i)=currSdItem;
        weiClus(i)=currWclus;
        weiItem(i)=currWitem;
        probTarg(i)=currProbTarg;
        llkAll(i)=llkAll(i-1);
    end    
    
    %% Determine trial params
    candAssign=currAssign;
    currTotalSd=sqrt((sdGlobe(i)^2)+(sdClus(i)^2)+(sdItem(i)^2));
    for j=1:numSubj
        candAssign2=candAssign(j,:);
        % Find likelihood of current assignments
        assProbs1=normpdf(guesses(candAssign(j,:),1,j),targCent(:,1,j)+(weiClus(i)*targClus(:,1,j))+(weiItem(i)*targItem(:,1,j)),currTotalSd);
        assProbs2=normpdf(guesses(candAssign(j,:),2,j),targCent(:,2,j)+(weiClus(i)*targClus(:,2,j))+(weiItem(i)*targItem(:,2,j)),currTotalSd);
        assProbs=assProbs1.*assProbs2;
        switch1=randsample(1:numI,1,true,((sum(assProbs)./(assProbs+(10^-300)))+(10^-300)));
        switch2=randsample(setdiff(1:numI,switch1),1,true,(sum(assProbs(1:numI~=switch1))./(assProbs(1:numI~=switch1)+(10^-300)))+(10^-300));
        sampSwitch2=fliplr([switch1 switch2]);
        candAssign2(1,[switch1 switch2])=candAssign2(1,sampSwitch2);
        %candAssign2=candAssign2(randperm(numI));
        llkAss=calcLikelihood4(targCent(:,:,j),targClus(:,:,j),targItem(:,:,j),guesses(:,:,j),sdGlobe(i),sdClus(i),sdItem(i),weiClus(i),weiItem(i),numC,corrClus,candAssign2,probTarg(i));
  
        a=prod(prod((llkAss+(10^-100))./(llkTotal(j,:)+(10^-100)),2));
        check=rand()<a;
        if a>1 || check
           assignments(j,:,i)=candAssign2; 
           llkTotal(j,:)=llkAss;
        else
           assignments(j,:,i)= candAssign(j,:);
        end
        
    end
    %assignments(:,:,i)=candAssign;
    
end

% sdGlobe(1:numBurn,:)=[];
% sdClus(1:numBurn,:)=[];
% sdItem(1:numBurn,:)=[];
% weiClus(1:numBurn,:)=[];
% weiItem(1:numBurn,:)=[];
% assignments(:,:,1:numBurn)=[];

    function input2=bound2(input,shift)
        input2=normrnd(input,shift);
        while input2<=0 || input2>=1
            input2=normrnd(input,shift);
        end
        
    end

    function input2=bound(input,shift)
        input2=normrnd(input,shift);
        if input2<=0 || input2>=1
            input2=input;
        end
    end

    function input2=bound3(input,shift)
        input2=normrnd(input,shift);
        if input2<=0 
            input2=input;
        end
    end
end







