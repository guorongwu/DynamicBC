function DynamicBC_ROI_dFC_mat2txt(matfile,save_dir)
% matfile  = 'F:\Data\TR645_AAL\FCM\s1_01\TV_s1_01_FCM.mat';
% save_dir = 'F:\Data\TR645_AAL\FCM\s1_01_txt\';
[fpath,name,~] = fileparts(matfile);
if nargin<2
    save_dir = fpath;
else
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
end
load(matfile)
nobs = length(FCM.Matrix);
num0 = ceil(log10(nobs))+2;
for i=1:nobs
    data = full(FCM.Matrix{i});
    num1 = num0 - length(num2str(i));
    fname = fullfile(save_dir,[name,repmat('0',1,num1),num2str(i),'.txt']);
    save(fname,'data','-ascii')
end