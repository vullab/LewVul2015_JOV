function pDist=findPairDist(currTarg)

% Cross Correlations
pDist=nan(size(currTarg,1),size(currTarg,1));
for m=1:size(currTarg,1)
    for n=1:size(currTarg,1)
        pd=sqrt(sum((currTarg(m,:)-currTarg(n,:)).^2));
        pDist(m,n)=pd;
    end
end



end