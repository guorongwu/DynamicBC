function DynamicBC_ROI_dBC_vis(matfile,fcec)
%% matfile = 'F:\Data\TR645_AAL\FCM\s1_01\TV_s1_01_FCM.mat'
%% fcec=1: FC, else: GC
load(matfile)
if fcec==1
    Matrix = FCM.Matrix;
else
    Matrix = GCM.Matrix;
end
n=length(Matrix);
for i=1:n
    f=Matrix{i};
    if i==1
        mn = min(f(:));
        mx = max(f(:));
    else
        mn = min([f(:);mn]);
        mx = max([f(:);mx]);
    end
end
figure(1);colormap(jet);

for i=1:n
    imagesc(Matrix{i});
    colorbar;%set(gca, 'CLim', [mn mx]);
    pause(0.15);
end
