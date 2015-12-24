function [params]=run_fitrelPos(targs,guesses,misAssign,clusAssign,fname)
% Use fminsearch to find the best fitting rho and theta parameters

if ~exist(fname)
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    
    rpSearch=@(x)relPosLlk(x,targs(1:30),guesses(:,1:30),misAssign(:,1:30),clusAssign(1:30));
    [xparam1 nll1]=fminsearchbnd(rpSearch,[.02 .5],[0 0],[5 pi],options);
    
    rpSearch=@(x)relPosLlk(x,targs(31:70),guesses(:,31:70),misAssign(:,31:70),clusAssign(31:70));
    [xparam2 nll2]=fminsearchbnd(rpSearch,[.02 .5],[0 0],[5 pi],options);
    
    params={xparam1 xparam2 ; nll1 nll2 };
    save(fname,'params');
else
    load(fname)
end
params=params(1,1:2);

    function logllk=relPosLlk(x,targs,guesses,misAssign,clusAssign)
        rho=x(1);
        theta=x(2);
        allLogllk=nan(size(targs));
        for ei=1:size(targs,2)
            % Extract clusters
            currTarg=targs{ei};
            currClus=clusAssign{ei};
            meanGloTarg=mean(currTarg);
            numC=length(unique(currClus));
            meanClusTarg=nan(numC,2);
            for ci=1:numC
                meanClusTarg(ci,:)=mean(currTarg(currClus==ci,:),1);
            end
            
            % Calculate relative distances
            radial=[];
            angular=[];
            cents=[repmat(meanGloTarg,numC,1);meanClusTarg(currClus,:)];
            projs=[meanClusTarg;currTarg];
            dif=cents-projs;
            [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
            radial=[radial; radTemp];
            angular=[angular; angTemp];
            
            for si=1:size(guesses,1)
                % Undo misassociations
                guesses2=guesses{si,ei};
                guesses3=guesses2(misAssign{si,ei},:);
                
                % Split into relative distances
                meanGloGuess=mean(guesses3);
                meanClusGuess=nan(numC,2);
                for ci=1:numC
                    meanClusGuess(ci,:)=mean(guesses3(currClus==ci,:),1);
                end
                cents=[repmat(meanGloGuess,numC,1);meanClusGuess(currClus,:)];
                projs=[meanClusGuess;guesses3];
                dif=cents-projs;
                [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
                
                radial2=[radial radTemp];
                angular2=[angular angTemp];
                
                rej=find(radial(:,1)==0);
                radial2(rej,:)=[];
                angular2(rej,:)=[];
                radial3=log10(radial2);
                rad_llk=normpdf(radial3(:,1),radial3(:,2),rho);
                ang_llk=circ_vmpdf(angular2(:,1),angular2(:,2),1/(theta^2));
                llk=rad_llk.*ang_llk+(10^-300);
                allLogllk(si,ei)=sum(log(llk));
            end
        end
        logllk=-sum(allLogllk(:));
    end

    function logllk=relPosLlk2(x,targs,guesses,misAssign,clusAssign)
        rho=x(1);
        theta=x(2);
        allLogllk=nan(size(targs));
        for ei=1:size(targs,2)
            % Extract clusters
            currTarg=targs{ei};
            currClus=clusAssign{ei};
            meanGloTarg=mean(currTarg);
            numC=length(unique(currClus));
            meanClusTarg=nan(numC,2);
            for ci=1:numC
                meanClusTarg(ci,:)=mean(currTarg(currClus==ci,:),1);
            end
            
            % Calculate relative distances
            radial=[];
            angular=[];
            cents=[repmat(meanGloTarg,numC,1);meanClusTarg(currClus,:)];
            projs=[meanClusTarg;currTarg];
            dif=cents-projs;
            [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
            radial=[radial; radTemp];
            angular=[angular; angTemp];
            
            for si=1:size(guesses,1)
                % Undo misassociations
                guesses2=guesses{si,ei};
                guesses3=guesses2(misAssign{si,ei},:);
                
                % Split into relative distances
                meanGloGuess=mean(guesses3);
                meanClusGuess=nan(numC,2);
                for ci=1:numC
                    meanClusGuess(ci,:)=mean(guesses3(currClus==ci,:),1);
                end
                cents=[repmat(meanGloGuess,numC,1);meanClusGuess(currClus,:)];
                projs=[meanClusGuess;guesses3];
                dif=cents-projs;
                [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
                
                radial2=[radial radTemp];
                angular2=[angular angTemp];
                
                rej=find(radial(:,1)==0);
                radial2(rej,:)=[];
                angular2(rej,:)=[];
                radial3=log10(radial2);
                rad_llk=normpdf(radial3(:,1),radial3(:,2),rho);
                ang_llk=circ_vmpdf(angular2(:,1),angular2(:,2),1/(theta^2));
                llk=rad_llk.*ang_llk+(10^-300);
                allLogllk(si,ei)=sum(log(llk));
            end
        end
        logllk=-mean(allLogllk(:));
    end

    function logllk=relPosLlk3(x,targs,guesses,misAssign,clusAssign)
        rho=x(1);
        theta=x(2);
        allLogllk=nan(size(targs));
        for ei=1:size(targs,2)
            % Extract clusters
            currTarg=targs{ei};
            currClus=clusAssign{ei};
            meanGloTarg=mean(currTarg);
            numC=length(unique(currClus));
            meanClusTarg=nan(numC,2);
            for ci=1:numC
                meanClusTarg(ci,:)=mean(currTarg(currClus==ci,:),1);
            end
            
            % Calculate relative distances
            radial=[];
            angular=[];
            cents=[repmat(meanGloTarg,numC,1);meanClusTarg(currClus,:)];
            projs=[meanClusTarg;currTarg];
            dif=cents-projs;
            [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
            radial=[radial; radTemp];
            angular=[angular; angTemp];
            
            for si=1:size(guesses,1)
                % Undo misassociations
                guesses2=guesses{si,ei};
                guesses3=guesses2(misAssign{si,ei},:);
                
                % Split into relative distances
                meanGloGuess=mean(guesses3);
                meanClusGuess=nan(numC,2);
                for ci=1:numC
                    meanClusGuess(ci,:)=mean(guesses3(currClus==ci,:),1);
                end
                cents=[repmat(meanGloGuess,numC,1);meanClusGuess(currClus,:)];
                projs=[meanClusGuess;guesses3];
                dif=cents-projs;
                [angTemp radTemp]=cart2pol(dif(:,1),dif(:,2));
                
                radial2=[radial radTemp];
                angular2=[angular angTemp];
                
                rej=find(radial(:,1)==0);
                radial2(rej,:)=[];
                angular2(rej,:)=[];
                radial3=log10(radial2);
                rad_llk=normpdf(radial3(:,1),radial3(:,2),rho);
                ang_llk=circ_vmpdf(angular2(:,1),angular2(:,2),1/(theta^2));
                llk=rad_llk.*ang_llk+(10^-300);
                allLogllk(si,ei)=sum(log(llk));
            end
        end
        logllk=-median(allLogllk(:));
    end
end





