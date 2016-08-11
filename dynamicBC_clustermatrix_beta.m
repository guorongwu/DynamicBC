function [IDX_subj,IDX_subjre]=dynamicBC_clustermatrix_beta(k,outputd,matname,subjdir,CluMet)
evalstr = ['mats = ',matname,';'];
dirpwd = pwd;
% subjdir = spm_select(1,'dir','Subj dir');
SubFold = dir(subjdir);
NumOfSubFold = size(SubFold,1)-2;
matrixall = [];
inums = 1;
for i = 1:NumOfSubFold
    if isdir([subjdir,SubFold(i+2).name])
        cd([subjdir,SubFold(i+2).name]);
        matsnames = dir('*.mat');
        load(matsnames.name)
        eval(evalstr);
        nmats = length(mats);
        for imat = 1:nmats
            matrixs(:,imat) = reshape(mats{imat},prod(size(mats{imat})),1);
        end
        matrixs = full(matrixs);
        matrixall = [matrixall,matrixs];
        matrix_all{inums} = matrixs;
        subjlist{inums} = SubFold(i+2).name;
        inums = inums+1;
    end
end
dims = size(mats{imat});
%%

%%
[IDX,C,sumd,D] = kmeans(matrixall',k,'distance',CluMet);
mkdir([outputd,'mat_',CluMet,'_Kmeans_',num2str(k)]);
for i = 1:k
    indtemp = find(IDX==i);
    CIND_group{i,1} = indtemp;
    meanmaptemp = mean(matrixall(:,indtemp),2);
    DAT = reshape(meanmaptemp,dims(1),dims(2));
    save([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,'Cluster_',num2str(i),'.mat'],'DAT')
    DAT0(:,i) = mean(matrixall(:,indtemp),2);
end

for isubj = 1:inums-1
    DATUSED = matrix_all{isubj};
    
    [IDX,C,sumd,D] = kmeans(DATUSED',k,'distance',CluMet);
    mkdir([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,subjlist{isubj}]);
    for i = 1:k
        indtemp = find(IDX==i);
%         CINDFLS{k,i,isubj} = indtemp;
        meanmaptemp = mean(DATUSED(:,indtemp),2);
        DAT = reshape(meanmaptemp,dims(1),dims(2));
        save([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,...
            subjlist{isubj},filesep,'Cluster',num2str(i),'.mat'],'DAT')
%         DAT = DAT.*MASKdat;
%         rest_WriteNiftiImage(DAT,Header,[outputd,'Corr_Kmeans_',num2str(k),...
%             filesep,SubFold(isubj+2).name,filesep,'Cluster',num2str(i),'.nii']);
    end
    IDX_subj(isubj,:) = IDX;
end


for isubj = 1:inums-1
    for i = 1:k
        dat1 = load([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),...
            filesep,subjlist{isubj},filesep,'Cluster',num2str(i),'.mat']);
        DAT1(:,i) = reshape(dat1.DAT,prod(dims),1);
    end
    DATrebuild{isubj,1} = DAT1;
end
% save test0001
for i = 1:inums-1
    perm_reb = perms(1:k);
    for iperm = 1:size(perm_reb)
        DATtemp = DATrebuild{i};
        
        [r p] = corr(DAT0,DATtemp(:,perm_reb(iperm,:)));
        R = diag(r);
        Rtemp(iperm) = sum(R);
    end
    maxR = find(Rtemp==max(Rtemp));
    maxperm{i} = perm_reb(maxR,:);
end

for i = 1:inums-1
    mkdir([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,'sorted',subjlist{i}]);
    maxpermtemp = maxperm{i};
    IDX0 = IDX_subj(i,:);
    IDX01 = zeros(size(IDX0));
    for j = 1:k
        copyfile([outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,subjlist{i},filesep,'Cluster',num2str(maxpermtemp(j)),'.mat'],...
            [outputd,'mat_',CluMet,'_Kmeans_',num2str(k),filesep,'sorted',subjlist{i},filesep,'Cluster',num2str(j),'_',num2str(maxpermtemp(j)),'.mat']);
        IDX01(IDX0==maxpermtemp(j)) = j;
    end
    IDX_subjre(i,:) = IDX01;
end
cd(dirpwd);

