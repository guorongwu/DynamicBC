function [CC,C,subid,str] = DynamicBC_export_mat(datadir,filtstr,id)
% filtstr='MedianCluster1';
% filtstr='FCM';
% filtstr='GCM';
% filtstr='CM';
 
% datadir = 'F:\Data\TR645_AAL\FCM';
% datadir = 'F:\Data\TR645_GC_AAL\GCM';
% datadir = 'F:\Data\TR645_GC_AAL\GCM_clu'; 

matfile0 = spm_select('FPListRec',datadir,[filtstr,'.*\.mat$']);
C=[];str={};
nfile = size(matfile0,1);
for i=nfile:-1:1
    a = strcat(matfile0(i,:));
    [~,subid{i,1}] = fileparts(fileparts(a)); 
    fprintf('%s\n',a)
    b = load(a);
    flag=isfield(b,{'FCM','GCM','DAT'});
    if flag(1)
        c = b.FCM.variance;
    elseif flag(2)
        c = b.GCM.variance;
    elseif flag(3)
        c = b.DAT;
    end
    
    if nargin<3%~exist('id','var')
        if i==nfile
            tt = c-c'; tt(tt<10^-6)=0;
            if any(tt(:))
                id = find([ones(size(tt))-eye(size(tt))]);
                disp('asymmetry matrix')
                flag2 = 1;
            else
                id = find(tril(ones(size(tt)),-1));
                disp('symmetry matrix')
                flag2 = 0;
            end
            k=1;
            for x=1:size(c,1)
                for y = 1:size(c,2)
                    if x<y
                        str{1,k} = [num2str(x),'_to_',num2str(y)]; k=k+1;
                    end
                    if x>y&&flag2
                        str{1,k} = [num2str(x),'_to_',num2str(y)]; k=k+1;
                    end
                end
            end
        end
        C(i,:)= c(id);
    else
        id2  = sub2ind(size(c),id(:,1),id(:,2));
        C(i,:)= c(id2);
        k=1;
        for j=1:size(id,1)
            str{1,k} = [num2str(id(j,1)),'_to_',num2str(id(j,2))]; k=k+1;
        end
    end
end
CC = [{'Subject'},str;subid num2cell(C)];