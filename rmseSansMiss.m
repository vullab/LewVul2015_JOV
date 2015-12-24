function [rmse2]=rmseSansMiss(targs,guesses,allAssign)
% Calculates RMSE without misassociations
rmse2=nan(size(targs));
for si=1:size(targs,1)
    for ei=1:size(targs,2)
        
        tdist=mean(sqrt(sum((targs{si,ei}-guesses{si,ei}(allAssign{si,ei},:)).^2,2)));
        rmse2(si,ei)=tdist;
    end
end


