function [params]=run_fitrelPos2(targs,guesses,misAssign,clusAssign,fname)
% Use fminsearch to find the best fitting rho and theta parameters
% 12/8/2014-V2-Fitting for each structure
if ~exist(fname)
    options = optimset('MaxFunEvals', 1000, 'MaxIter', 1000);
    for ei1=1:7
        inds=(1+(10*(ei1-1))):10*ei1;
        rpSearch=@(x)relPosLlk(x,targs(inds),guesses(:,inds),misAssign(:,inds),clusAssign(inds));
        [xparam1 nll1]=fminsearchbnd(rpSearch,[.02 .5],[0 0],[5 pi],options);
        
        params{ei1}={xparam1  nll1  };
    end
    save(fname,'params');
else
    load(fname)
end

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





