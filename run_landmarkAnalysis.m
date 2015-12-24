%% Landmark analysis

numSubj=size(guesses,1);
clusErrX=[];
clusErrY=[];
clusPosX=[];
clusPosY=[];
clusType=[];
for si=1:numSubj
    count=1;
    for ei=1:70
        currC=numClu(ceil(ei/10));
        currIn=numIns(ceil(ei/10));
        
        currT=targs{si,ei};
        guesses2=guesses{si,ei}(subjAssign{si,ei},:);
        currG=guesses2;
        
        currErr=abs((currT)-(currG));
        
        clusErrX(count:count+(size(currT,1)-1),si)=(currErr(:,1));
        clusErrY(count:count+(size(currT,1)-1),si)=(currErr(:,2));
        
        
        clusPosX(count:count+(size(currT,1)-1),si)=currT(:,1);
        clusPosY(count:count+(size(currT,1)-1),si)=currT(:,2);
        
        clusType(count:count+(size(currT,1)-1),1)=ceil(ei/10);
        count=count+size(currT,1);
    end
end