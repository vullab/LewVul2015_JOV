function llk=calcLikelihood4(targCent,targClus,targItem,guess,sdGlobal,sdClus,sdItem,weiClus,weiItem,numC,corrClus,assign,probTarg)

numI=size(targItem,1);
numIn=numI/numC;
numSubj=size(targItem,3);
% Calculate covariance matrix
varGlobal=sdGlobal^2;
varClus=sdClus^2;
varItem=sdItem^2;
prime=[2 3 5 7 11 13 17 19];
% Pdf for x
llk=nan(numSubj,3);
for i=1:numSubj
    covar=repmat(varGlobal,numI)+eye(numI)*varItem;
    covar(mod(sqrt((prime(corrClus(assign(i,:)))'*prime(corrClus(assign(i,:))))),1)==0)=covar(mod(sqrt((prime(corrClus(assign(i,:)))'*prime(corrClus(assign(i,:))))),1)==0)+varClus;
    llk(i,1)=mvnpdf(guess(assign(i,:),1,i),(targItem(:,1,i)*weiItem)+(targClus(:,1,i)*weiClus)+targCent(:,1,i),covar);
    llk(i,2)=mvnpdf(guess(assign(i,:),2,i),(targItem(:,2,i)*weiItem)+(targClus(:,2,i)*weiClus)+targCent(:,2,i),covar);
    
    targNorm=(1-probTarg)/(numI-1);
    assProbs=ones(numI,1);
    assProbs2=prod(assProbs,2);
    
    assProbs2(~(assign(i,:)==1:numI))=(targNorm).*assProbs2(~(assign(i,:)==1:numI),1);
    assProbs2((assign(i,:)==1:numI))=probTarg*assProbs2((assign(i,:)==1:numI),1);

    llk(i,3)=prod(assProbs2(:));
end
end