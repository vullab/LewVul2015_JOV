function totalDisp=dispEnv(targs,numClu)

totalDisp=nan(size(targs,1),size(targs,2),2);
for i=1:size(targs,1)
    for j=1:size(targs,2)
        totalDisp(i,j,1)=mean(std(targs{i,j},0,1),2);
        stdx=std(squeeze(mean(reshape(targs{i,j}(:,1)',size(targs{i,j},1)/numClu(ceil(j/10)),1,numClu(ceil(j/10))),1)),0,1);
        stdy=std(squeeze(mean(reshape(targs{i,j}(:,2)',size(targs{i,j},1)/numClu(ceil(j/10)),1,numClu(ceil(j/10))),1)),0,1);
        totalDisp(i,j,2)=mean([stdx stdy],2);
    end
end

end