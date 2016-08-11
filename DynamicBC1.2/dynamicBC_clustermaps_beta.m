function [IDX_subj,IDX_subjre]=dynamicBC_clustermaps_beta(k,outputd,Maskimg,subjdir,CluMet)
% answers = inputdlg('Cluster number','Cluster number',1,{'6'});
% k = str2num(answers{1});clear answers;
% outputd = spm_select(1,'dir','Outputdir');
% Maskimg = spm_select(1,'image','mask image');
MASKdat = rest_ReadNiftiImage(Maskimg);
maskind = find(MASKdat);
dims = size(MASKdat);
% subjdir = spm_select(1,'dir','Subj dir');
SubFold = dir(subjdir);
NumOfSubFold = size(SubFold,1)-2;
DATAorig = [];
DATAused = [];
for i = 1:NumOfSubFold
    dirtemp = [subjdir,SubFold(i+2).name];
    [AllV,VoxS,ImgFL,Header,nVolumn] = rest_to4d(dirtemp);
    Dat = reshape(AllV,prod(dims),nVolumn);
    DATAorig = [DATAorig,Dat];
    DATAused = [DATAused,Dat(maskind,:)];
%     Nlength = nVolumn;
    DATAorigSubj{i} = Dat;
    DATAusedSubj{i} = Dat(maskind,:);
end
%%
[IDX,C,sumd,D] = kmeans(DATAused',k,'distance',CluMet);
mkdir([outputd,CluMet,'_Kmeans_',num2str(k)]);
for i = 1:k
    indtemp = find(IDX==i);
    CIND_group{i,1} = indtemp;
    meanmaptemp = mean(DATAorig(:,indtemp),2);
    DAT = reshape(meanmaptemp,dims(1),dims(2),dims(3));
    DAT = DAT.*MASKdat;
    rest_WriteNiftiImage(DAT,Header,[outputd,CluMet,'_Kmeans_',num2str(k),filesep,'Cluster_',num2str(i),'.nii']);
    DAT0(:,i) = mean(DATAused(:,indtemp),2);
end

for isubj = 1:NumOfSubFold
%     DATAorigSubj{i} = DAT;
%     DATAusedSubj{i} = Dat(maskind,:);
    DATORIG = DATAorigSubj{isubj};
    DATUSED = DATAusedSubj{isubj};
    
    [IDX,C,sumd,D] = kmeans(DATUSED',k,'distance',CluMet);
    mkdir([outputd,CluMet,'_Kmeans_',num2str(k),filesep,SubFold(isubj+2).name]);
    for i = 1:k
        indtemp = find(IDX==i);
        CINDFLS{k,i,isubj} = indtemp;
        meanmaptemp = mean(DATORIG(:,indtemp),2);
        DAT = reshape(meanmaptemp,dims(1),dims(2),dims(3));
        DAT = DAT.*MASKdat;
        rest_WriteNiftiImage(DAT,Header,[outputd,CluMet,'_Kmeans_',num2str(k),...
            filesep,SubFold(isubj+2).name,filesep,'Cluster',num2str(i),'.nii']);
    end
    IDX_subj(isubj,:) = IDX;
end

% for i = 1:6
%     dat = rest_ReadNiftiImage(['D:\DBC\TestForGroup\Corr_FLS_Kmeans_6\Test',num2str(i),'.nii']);
%     DAT0(:,i) = dat(maskind);
% end
for isubj = 1:NumOfSubFold
    for i = 1:k
        dat1 = rest_ReadNiftiImage([outputd,CluMet,'_Kmeans_',num2str(k),...
            filesep,SubFold(isubj+2).name,filesep,'Cluster',num2str(i),'.nii']);
        DAT1(:,i) = dat1(maskind);
    end
    DATrebuild{isubj,1} = DAT1;
end

for i = 1:NumOfSubFold
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

for i = 1:NumOfSubFold
    mkdir([outputd,CluMet,'_Kmeans_',num2str(k),filesep,'sorted',SubFold(i+2).name]);
    maxpermtemp = maxperm{i};
    IDX0 = IDX_subj(i,:);
    IDX01 = zeros(size(IDX0));
    for j = 1:k
        copyfile([outputd,CluMet,'_Kmeans_',num2str(k),filesep,SubFold(i+2).name,filesep,'Cluster',num2str(maxpermtemp(j)),'.nii'],...
            [outputd,CluMet,'_Kmeans_',num2str(k),filesep,'sorted',SubFold(i+2).name,filesep,'Cluster',num2str(j),'_',num2str(maxpermtemp(j)),'.nii']);
        IDX01(IDX0==maxpermtemp(j)) = j;
    end
    IDX_subjre(i,:) = IDX01;
end
    