clc,clear
% if 1
    matname ={
'F:\Data\TR645_AAL\FCM_clu\mat_sqeuclidean_Kmeans_4\Centroid_1.mat'
'F:\Data\TR645_AAL\FCM_clu\mat_sqeuclidean_Kmeans_4\Centroid_2.mat'
'F:\Data\TR645_AAL\FCM_clu\mat_sqeuclidean_Kmeans_4\Centroid_3.mat'
'F:\Data\TR645_AAL\FCM_clu\mat_sqeuclidean_Kmeans_4\Centroid_4.mat'
    };


n=length(matname);
for i=1:n
    clear DAT
    load(matname{i})
    if i==1
        mn = min(DAT(:));
        mx = max(DAT(:));
    else
        mn = min([DAT(:);mn]);
        mx = max([DAT(:);mx]);
    end
end
scrsz = get(groot,'ScreenSize');
figure('Position',[100 scrsz(4)/3 scrsz(3)*4/5 scrsz(4)/4])
for i=1:n
    clear DAT
    load(matname{i})
    [fpath,name]=fileparts(matname{i});
    subplot(1,n,i);imagesc(DAT);
    axis square;axis off
    set(gca, 'CLim', [mn mx]);
    pause(0.15);
    title(strrep(name,'_','-'))
end
colormap(jet);
