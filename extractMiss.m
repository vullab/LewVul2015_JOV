function [pClus pItem]=extractMiss(assignments,corrClus)

numSamps=size(assignments,3);
clusters=unique(corrClus);
pClus=nan(numSamps,1);
pItem=nan(numSamps,1);
corrClus2=repmat(corrClus,size(assignments,1),1);
corrOrder=repmat(1:size(assignments,2),size(assignments,1),1);
for i=1:numSamps
    currAssign=assignments(:,:,i);
    misMatch=corrClus(currAssign(:,:))~=corrClus2;
    pClus(i)=sum(misMatch(:))/numel(currAssign);
    
    switchCond=0;
    for j=1:length(clusters)
        inds=(corrClus2==clusters(j)).*(~misMatch);
        switchCond=switchCond+(sum(corrOrder(logical(inds))~=currAssign(logical(inds)))/length(corrOrder(logical(inds))));
    end
    switchCond=switchCond/length(clusters);
    pItem(i)=switchCond;
end












