function [] = DynamicBC_fls_FC(ROI_sig,mu,save_info)
[nobs,nvar] = size(ROI_sig);
ROI_sig = zscore(ROI_sig);
num0 = ceil(log10(nobs))+2;
% mu=100;
% tic
nii_name = cell(nobs,1);
if save_info.flag_nii % save nii: seed FC,FCD.
    data_save = zeros(save_info.v.dim);
    v = save_info.v;
    v.fname =  strcat(save_info.nii_mat_name);
    [pathstr, name, ext] = fileparts(v.fname) ;
    
    if ~save_info.flag_1n % FCD    beta = zeros(nvar,nobs); % FCD
        pathstr_fcd = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FCD_FLS_map',filesep));
        try
            mkdir(pathstr_fcd)
        end
        beta = zeros(nvar,nobs);
        parfor k=1:nvar
            betak = zeros(nvar,nobs);
            for j=1:nvar
                if j~=k
                    betak(j,:) = wgr_fls(ROI_sig(:,k), ROI_sig(:,j), mu);
                end
            end
            beta(k,:) = DynamicBC_nansum(abs(betak),1); 
        end
        
        for k=1:nobs
            num1 = num0 - length(num2str(k));
            v.fname = fullfile(pathstr_fcd,[name,repmat('0',1,num1),num2str(k),ext]);
            data_save(save_info.index) = beta(:,k); 
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;
        end
    else %if save_info.flag_1n % seed FC
        pathstr_fc = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FC_FLS_map',filesep));
        try
            mkdir(pathstr_fc)
        end
        Matrix = zeros(nvar,nobs);
        save_info.seed_signal = zscore(save_info.seed_signal);
        for k=1:nvar           
            Matrix(k,:) = wgr_fls(save_info.seed_signal, ROI_sig(:,k), mu);
        end
        for k=1:nobs
            num1 = num0 - length(num2str(k));
            v.fname = fullfile(pathstr_fc,[name,repmat('0',1,num1),num2str(k),ext]);
            data_save(save_info.index) = Matrix(:,k); 
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;
        end
    end
    
else % save as mat file
    
    beta = zeros(nvar,nvar,nobs); % FC
    for k=1:nvar
        betak = zeros(nvar,nobs);
        for j=1:nvar
%             fprintf('.')
            if j~=k
                betak(j,:) = wgr_fls(ROI_sig(:,k), ROI_sig(:,j), mu);
            end
        end
        beta(k,:,:) = betak;
        fprintf('.')
    end
    fprintf('.\n')
    for k=1:nobs
        FCM.Matrix{k} = (beta(:,:,k)+beta(:,:,k)')/2;
    end
%     save(save_info.nii_mat_name,'FCM');
end
% toc

%% variance calculation

if save_info.flag_nii % save nii: seed FC,FCD.
    str_typ = '_variance';
    if ~save_info.flag_1n % FCD
        pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FCD_FLS_Variance',filesep));
        try
            mkdir(pathstr_v)
        end
        data0 = zeros(nobs, nvar);
        for k=1:nobs
            v = spm_vol(nii_name{k,1}.fname);
            dat = spm_read_vols(v);
            data0(k,:) = dat(save_info.index);
        end
        data = var(data0,0,1);
        v.fname = fullfile(pathstr_v,[name,'_variance',ext]);
        data_save(save_info.index) = data; 
        spm_write_vol(v,data_save);
        
        datas = std(data0,0,1);
        v.fname = fullfile(pathstr_v,[name,'_std',ext]);
        data_save(save_info.index) = datas; 
        spm_write_vol(v,data_save);
        
        datam = mean(data0,1);
        v.fname = fullfile(pathstr_v,[name,'_mean',ext]);
        data_save(save_info.index) = datam; 
        spm_write_vol(v,data_save);
        
        cv = datas./datam;
        v.fname = fullfile(pathstr_v,[name,'_CV',ext]);
        data_save(save_info.index) = cv; 
        spm_write_vol(v,data_save);
        
        v.fname = fullfile(pathstr_v,[name,'_CV_abs',ext]);
        data_save(save_info.index) = abs(cv); 
        spm_write_vol(v,data_save);
        
    else
        pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FC_FLS_Variance',filesep));
        try
            mkdir(pathstr_v)
        end
        data0 = zeros(nobs, nvar);
        for k=1:nobs
            v = spm_vol(nii_name{k,1}.fname);
            dat = spm_read_vols(v);
            data0(k,:) = dat(save_info.index);
        end
        
        data = var(data0,0,1);
        v.fname = fullfile(pathstr_v,[name,'_variance',ext]);
        data_save(save_info.index) = data; 
        spm_write_vol(v,data_save);
        
        datas = std(data0,0,1);
        v.fname = fullfile(pathstr_v,[name,'_std',ext]);
        data_save(save_info.index) = datas; 
        spm_write_vol(v,data_save);
        
        datam = mean(data0,1);
        v.fname = fullfile(pathstr_v,[name,'_mean',ext]);
        data_save(save_info.index) = datam; 
        spm_write_vol(v,data_save);
        
        cv = datas./datam;
        v.fname = fullfile(pathstr_v,[name,'_CV',ext]);
        data_save(save_info.index) = cv; 
        spm_write_vol(v,data_save);
        
        v.fname = fullfile(pathstr_v,[name,'_CV_abs',ext]);
        data_save(save_info.index) = abs(cv); 
        spm_write_vol(v,data_save);
    end
    
else
    
    data_var = zeros(nobs, nvar*nvar);
    for k=1:nobs
        tmp_matrix  = FCM.Matrix{k};
        data_var(k,:) = tmp_matrix(:);               
    end
    data_var = var(data_var,0,1);
    data_var = reshape(data_var,nvar,nvar);
    FCM.variance = data_var;
    
    save_info.nii_mat_name = strrep(save_info.nii_mat_name,save_info.save_dir,fullfile(save_info.save_dir,'FCM',filesep));
    [fcm_dir] = fileparts(save_info.nii_mat_name);
    try
        mkdir(fcm_dir)
    end
    save(save_info.nii_mat_name,'FCM');
end

function b_FLS = wgr_fls(X,y,mu);
[N,K] = size(X);
G=zeros(N*K,N);
A=zeros(N*K,N*K);
mui = mu*eye(K);
ind = 1:K;

for i=1:N
    G(ind,i) = X(i,:);
    if i==1
        Ai = X(i,:)'*X(i,:) + mui;
        A(ind,ind)= Ai ;
        A(ind,ind+K)= - mui;
    elseif i~=1 && i~=N
        Ai = X(i,:)'*X(i,:) + 2*mui;
        A(ind,ind)= Ai ;
        A(ind,ind+K)= - mui;
        A(ind,ind-K)= - mui;
    else% i==N
        Ai = X(i,:)'*X(i,:) + mui;
        A(ind,ind)= Ai ;
        A(ind,ind-K)= - mui;
    end        
    ind = ind+K;
end
b_FLS = linsolve(A,G*y);
b_FLS = reshape(b_FLS,K,N)';

function y = DynamicBC_nansum(x,dim)
%NANSUM Sum, ignoring NaNs.
%   Y = NANSUM(X) returns the sum of X, treating NaNs as missing values.
%   For vector input, Y is the sum of the non-NaN elements in X.  For
%   matrix input, Y is a row vector containing the sum of non-NaN elements
%   in each column.  For N-D arrays, NANSUM operates along the first
%   non-singleton dimension.
%
%   Y = NANSUM(X,DIM) takes the sum along dimension DIM of X.
%
%   See also SUM, NANMEAN, NANVAR, NANSTD, NANMIN, NANMAX, NANMEDIAN.

%   Copyright 1993-2004 The MathWorks, Inc.


% Find NaNs and set them to zero.  Then sum up non-NaNs.  Cols of all NaNs
% will return zero.
x(isnan(x)) = 0;
if nargin == 1 % let sum figure out which dimension to work along
    y = sum(x);
else           % work along the explicitly given dimension
    y = sum(x,dim);
end
