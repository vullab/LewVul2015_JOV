%% run_pairDist_CrossCorr
% Put pairDist data and cross corr data into compatible forms
pairData=[];
ccData=[];
sameClus=[];

numC=[4 2 1 8 4 2 1];
numI=[4 4 4 8 8 8 8]./numC;
numI2=[4 4 4 8 8 8 8];
for si=1:35
    count=1;
    for ei=1:70
        
        currPair=subjPairDist{si,ei};
        currCross=subjCC2{si,ei};
        
        ei2=ceil(ei/10);
        sameC=zeros(numI2(ei2));
        for ci=1:numC(ei2)
            rr=(ci-1)*numI(ei2)+1:(ci-1)*numI(ei2)+numI(ei2);
            sameC(rr,rr)=1;
        end

        
        for m=1:size(currPair,1)
            for n=1:size(currPair,1)
                if m>n
                    sameClus(si,count)=sameC(m,n);
                    pairData(si,count)=currPair(m,n);
                    ccData(si,count)=currCross(m,n);
                    count=count+1;
                end
            end
        end
        
    end
end

% figure;plot(mean(pairData),mean(ccData),'k.');xlabel('Distance');ylabel('Error similarity')
% 
% figure;
% pairDataClus=pairData(:,logical(mean(sameClus)>0));
% pairDataNoClus=pairData(:,logical(mean(sameClus)==0));
% ccDataClus=ccData(:,logical(mean(sameClus)>0));
% ccDataNoClus=ccData(:,logical(mean(sameClus)==0));
% hold on
% plot(mean(pairDataClus),mean(ccDataClus),'r.');
% plot(mean(pairDataNoClus),mean(ccDataNoClus),'b.');
% xlabel('Distance');ylabel('Error similarity')
% hold off
% 
% 
% [a b]=regress(mean(ccDataClus)',[mean(pairDataClus)' ones(length(mean(pairDataClus)'),1)]);
% 
% [a2 b2]=regress(mean(ccDataNoClus)',[mean(pairDataNoClus)' ones(length(mean(pairDataNoClus)'),1)]);
% 
% 



