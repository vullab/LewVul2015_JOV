%% run_crossCorr

if ~exist(fullfile(strcat('subjCC3','.mat')))
    subjCC2=cell(length(subjid),size(distances,2));
    for i=1:length(subjid)
        for j=1:size(distances,2)
            currC=numClu(ceil(j/numEnv));
            currGuess=guesses{i,j};
            currTarg=targs{i,j};
            cCorr=findCrossCorr3(currGuess,currTarg,currC);
            subjCC2{i,j}=cCorr;
            
        end
    end
    save(fullfile(strcat('subjCC3','.mat')),'subjCC2')
else
    load(fullfile(strcat('subjCC3','.mat')))
end

% Get random distribution
randresp=rand(1000,2,2);
randresp=randresp-.5;
randCcorr=nan(1000,1);
for i=1:1000
    randCcorr(i)=(randresp(i,:,1)*randresp(i,:,2)')/(norm(randresp(i,:,1))*norm(randresp(i,:,2)));
end

sampling=.71;

subjClus=nan(size(subjCC2,3),7);
subjClus2=nan(size(subjCC2,3),7);
numC=[4 2 1 8 4 2 1];
numI=[4 4 4 8 8 8 8]./numC;
for si=1:size(subjCC2,1)
    for ei=1:7
        numO=size(subjCC2{1,(ei-1)*10+1},1);
        sel=zeros(numO,numO);
        for ci=1:numC(ei)
            rr=(ci-1)*numI(ei)+1:(ci-1)*numI(ei)+numI(ei);
            sel(rr,rr)=1;
        end
        
        sel(logical(eye(numO)))=0;
        temp=nan(10,1);
        temp2=nan(10,1);
        for ei2=1:10
            
            temp(ei2)=mean(subjCC2{si,(ei-1)*10+ei2}(logical(sel)));
            temp2(ei2)=mean(subjCC2{si,(ei-1)*10+ei2}(logical(~sel)));
        end
        
        subjClus(si,ei)=mean(temp);
        subjClus2(si,ei)=mean(temp2);
    end
end