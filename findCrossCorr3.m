function cCorr=findCrossCorr3(currGuess,currTarg,numC)
% v2-Modified to work with any clustering
numIn=size(currTarg,1)/numC;
currDifs=currGuess-currTarg;

% Cross Correlations
cCorr=nan(size(currGuess,1),size(currGuess,1));
for m=1:size(currGuess,1)
    for n=1:size(currGuess,1)
        total=(currDifs(m,:)*currDifs(n,:)')/(norm(currDifs(m,:))*norm(currDifs(n,:)));
        cCorr(m,n)=total;
    end
end



end