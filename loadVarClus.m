%% Load and process varClus data
fpath=fullfile(cd,'data');
load(fullfile(fpath,'GenEnvData10env3Used.mat'))
fname='varClusData'; % varClusData.txt-AKA data6_3
datafile=fopen(fullfile(fpath,strcat(fname,'.txt')));
data =textscan(datafile,'%s %f %f %f %f %f %f %s %f %f %f', 'delimiter',';');

taskid=data{2};

% Select experiment iteration
idcheck=taskid==6.3;

id=data{1}(idcheck);
count=data{3}(idcheck);
phase=data{10}(idcheck);
phase=phase+1;
xsel=data{4}(idcheck);
ysel=data{5}(idcheck);
xtarg=data{6}(idcheck);
ytarg=data{7}(idcheck);
item=data{8}(idcheck);
time=data{9}(idcheck);
setId=data{11}(idcheck);

subjid=unique(id);
remove={'E104BW' 'D126AJ','C126CK','C130YT','test1','testTestTEST'}; % 'D135WL','D131JK','D130VV'
% E104BW-Encountered an error, object did not appear
% D126AJ-Used tapping mnemonic
subjid=setdiff(subjid,remove);
phases=unique(phase);

% Actual positions
if ~exist(fullfile('procData.mat'))
    distances=nan(length(subjid),length(phases));
    centerDists=nan(length(subjid),length(phases));
    guesses=cell(length(subjid),length(phases));
    targs=cell(length(subjid),length(phases));
    relTargs=cell(length(subjid),length(phases));
    clusCenters=cell(length(subjid),length(phases));
    for i=1:length(subjid)
        for j=1:length(phases)
            
            currids=logical(strcmp(id,subjid{i}).*(phase==phases(j)));
            items=item(currids);
            uniItems=setdiff(unique(items),'[');
            xselcurr=xsel(currids);
            yselcurr=ysel(currids);
            xtargcurr=xtarg(currids);
            ytargcurr=ytarg(currids);
            
            dists=nan(length(uniItems),1);
            cdists=nan(length(uniItems),1);
            xyPos=nan(length(uniItems),2);
            orderTargs=roundsd(allStore{j,1}{1},6)';
            invCount=0;
            xyGuess=nan(size(orderTargs,1),2);
            xyTarg=nan(size(orderTargs,1),2);
            xyRel=nan(size(orderTargs,1),2);
            for k=1:length(uniItems)
                [mInd m]=find(strcmp(uniItems(k),items));
                currInd=mInd(end);
                match=sqrt(sum((orderTargs-repmat(roundsd([xtargcurr(currInd) ytargcurr(currInd)],6),size(orderTargs,1),1)).^2,2));
                findInd=match<1;
                
                xyGuess(findInd,:)=[xselcurr(currInd) yselcurr(currInd)];
                xyTarg(findInd,:)=[xtargcurr(currInd) ytargcurr(currInd)];
                xyRel(findInd,:)=[xtargcurr(currInd) ytargcurr(currInd)]-allStore{j,2}{1}(:,ceil(find(findInd)/allStore{j,4}))';
                
                dists(findInd)=sqrt(((xselcurr(currInd)-xtargcurr(currInd))^2)+((yselcurr(currInd)-ytargcurr(currInd))^2));
                cdists(findInd)=(((xselcurr(currInd)-xcent)^2)+((yselcurr(currInd)-ycent)^2));
            end
            clusRel=allStore{j,2}{1};
            
            guesses{i,j}=xyGuess;
            targs{i,j}=xyTarg;
            relTargs{i,j}=xyRel;
            clusCenters{i,j}=clusRel;
            centerDists(i,j)=mean(cdists);
            distances(i,j)=mean(dists);
            
        end
    end
    save(fullfile('procData.mat'),'guesses','targs','clusCenters','relTargs','subjid','numClu','numIns','numItems','numEnv','distances')
else
    load(fullfile('procData.mat'))
end
 numSets=size(targs,1);