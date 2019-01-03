function [varargout] = DynamicBC_sliding_window_GC(data,window,overlap,pvalue,order,save_info)
%% calculate sliding window functional connectivity (bivariate)
% overlap = 0.1; %e.g. time bin: [1:50],[46:95], [91:140],...
% window = 50;
[nobs, nvar] = size(data);
step=ceil(window-overlap*window); % 10% overlap
if ~step||step<0
    error('you must reset overlap size!');
end
if window>nobs
   fprintf('There are only %d time points < window size = %d',nobs, window)
   error('you must reset window size!');
end
if window==nobs
    slides=1;
else
%     slides=floor((nobs-window)/step);
    slides=floor((nobs-window)/step)+1;
end
num0 = ceil(log10(slides))+2;
GC = cell(slides,1); 
t1=1-step;
t2=window-step;
%sliding window
nii_name = cell(slides,4);
save_info.flag_par = 1; % default open par
for k=1:slides
    t1=t1+step;
    t2=t2+step;
%     if k==slides
%         t2 = nobs;
%     end
    disp([t1 t2])
    dat = data(t1:t2,:);
    if save_info.slw_alignment==1
        k1 = t1;
    else
        k1 = floor((t1+t2)./2);
    end
    if save_info.flag_nii % save nii: seed GC,GCD.
        v = save_info.v;
        v.fname = strcat(save_info.nii_mat_name);
        [pathstr, name, ext] = fileparts(v.fname) ;
        
        num1 = num0 - length(num2str(k1));
        data_save = zeros(save_info.v.dim);
        if ~save_info.flag_1n % GCD
            pathstr_gcd = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'GCD_map',filesep));
            if k==1
                try
                    mkdir(pathstr_gcd)
                end
            end
            flag_gcd=1;
            [GCM] = DynamicBC_PWGC_SLW(dat,nvar,save_info.flag_1n,flag_gcd,pvalue,order,[],save_info.flag_par);
            Matrix_bin = sparse(nvar,nvar);
            Matrix_wei = Matrix_bin;
            for i=1:GCM.nbin
                for j=1:GCM.nbin
                    Matrix_wei(GCM.indX{i},GCM.indX{j}) =  GCM.matrix{i,j};
                    Matrix_bin(GCM.indX{i},GCM.indX{j}) =  single(GCM.matrix{i,j}>=GCM.gc_th);
                end
            end
            Matrix_wei(1:nvar+1:end)=0; %remove diag value.
            Matrix_wei(isnan(Matrix_wei)) = 0;
            Matrix_bin(1:nvar+1:end)=0; %remove diag value.

            data_save(save_info.index) = sum(abs(Matrix_bin),1);
            v.fname = fullfile(pathstr_gcd,[name,'_IN_bin_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;
            
            data_save(save_info.index) = sum(abs(Matrix_wei),1);
            v.fname = fullfile(pathstr_gcd,[name,'_IN_wei_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,2} = v;

            data_save(save_info.index) = sum(abs(Matrix_bin),2);
            v.fname = fullfile(pathstr_gcd,[name,'_OUT_bin_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,3} = v;
            
            data_save(save_info.index) = sum(abs(Matrix_wei),2);
            v.fname = fullfile(pathstr_gcd,[name,'_OUT_wei_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,4} = v;
            
        else  %if save_info.flag_1n % seed GC
            pathstr_g = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'seed_GCmap',filesep));
            if k==1
                try
                    mkdir(pathstr_g)
                end
            end
            flag_gcd=0;
            
            [GCM] = DynamicBC_PWGC_SLW(dat,nvar,save_info.flag_1n, flag_gcd, pvalue,order,save_info.seed_signal(t1:t2,:),save_info.flag_par);
            
            Matrix = zeros(nvar,2);
            for i=1:GCM.nbin
                Matrix(GCM.indX{i},:) =  GCM.matrix{i,1};
            end
            data_save(save_info.index) = Matrix(:,1);
            v.fname = fullfile(pathstr_g,[name,'_OUT_',repmat('0',1,num1),num2str(k1),ext]);            
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;

            data_save(save_info.index) = Matrix(:,2);
            v.fname = fullfile(pathstr_g,[name,'_IN_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,2} = v;
            
        end

    else % ROIwise save as mat file
        flag_gcd=0;
        [GCM] = DynamicBC_PWGC_SLW(dat,nvar,save_info.flag_1n,flag_gcd, pvalue,order,[],save_info.flag_par);
        Matrix = zeros(nvar,nvar);
        for i=1:GCM.nbin
            for j=1:GCM.nbin
                Matrix(GCM.indX{i},GCM.indX{j}) =  GCM.matrix{i,j};
            end
        end
        Matrix(1:nvar+1:end)=0;
        GC{k,1}  = Matrix; 
        time_alignment(k) = k1; 
    end
end

%% variance calculation

if save_info.flag_nii % save nii: seed GC,GCD.
    if ~save_info.flag_1n % GCD
        pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'GCD_Variance',filesep));
        try
            mkdir(pathstr_v)
        end
        str_typ = {'_In_bin_variance','_In_wei_variance','_Out_bin_variance','_Out_wei_variance'};
        str_typ1 = {'_In_bin_std','_In_wei_std','_Out_bin_std','_Out_wei_std'};
        str_typ2 = {'_In_bin_mean','_In_wei_mean','_Out_bin_mean','_Out_wei_mean'};
        str_typ3 = {'_In_bin_CV','_In_wei_CV','_Out_bin_CV','_Out_wei_CV'};
        for j=1:4
            data0 = zeros(slides, nvar);
            for k=1:slides
                v = spm_vol(nii_name{k,j}.fname);
                dat = spm_read_vols(v);
                data0(k,:) = dat(save_info.index);
            end
            data = var(data0,0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ{j},ext]);
            data_save(save_info.index) = data; 
            spm_write_vol(v,data_save);
            
            datas = std(data0,0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ1{j},ext]);
            data_save(save_info.index) = datas; 
            spm_write_vol(v,data_save);
            
            datam = mean(data0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ2{j},ext]);
            data_save(save_info.index) = datam; 
            spm_write_vol(v,data_save);
            
            cv = datas./datam;
            v.fname = fullfile(pathstr_v,[name,str_typ3{j},ext]);
            data_save(save_info.index) = cv; 
            spm_write_vol(v,data_save);
        end
    else
        if slides>1 
            pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'GC_Variance',filesep));
            try
                mkdir(pathstr_v)
            end
            for j=1:2
                data0 = zeros(slides, nvar);
                for k=1:slides
                    v = spm_vol(nii_name{k,j}.fname);
                    dat = spm_read_vols(v);
                    data0(k,:) = dat(save_info.index);
                end
                
                if j==1
                    postfix = '_Out';
                else
                    postfix = '_IN';
                end
                data = var(data0,0,1);
                v.fname = fullfile(pathstr_v,[name,postfix,'_variance',ext]);
                data_save(save_info.index) = data; 
                spm_write_vol(v,data_save);
                
                datas = std(data0,0,1);
                v.fname = fullfile(pathstr_v,[name,postfix,'_std',ext]);
                data_save(save_info.index) = datas; 
                spm_write_vol(v,data_save);
                
                datam = mean(data0,1);
                v.fname = fullfile(pathstr_v,[name,postfix,'_mean',ext]);
                data_save(save_info.index) = datam; 
                spm_write_vol(v,data_save);
                
                v.fname = fullfile(pathstr_v,[name,postfix,'_CV',ext]);
                data_save(save_info.index) = datas./datam; 
                spm_write_vol(v,data_save);
            end
        end
    end
    
else
    
    data_var = zeros(slides, nvar*nvar);
    for k=1:slides
        tmp_matrix  = GC{k,1};
        data_var(k,:) = tmp_matrix(:); 
    end
    data_var = var(data_var,0,1);
    data_var = reshape(data_var,nvar,nvar);
    
end

if ~save_info.flag_nii
    pvalue = GCM.pvalue;
    gc_th = GCM.gc_th;
    clear GCM
    GCM.pvalue = pvalue;
    GCM.gc_th = gc_th;
    GCM.Matrix = GC;
    GCM.variance = data_var;
    GCM.time_alignment = time_alignment;
    save_info.nii_mat_name = strrep(save_info.nii_mat_name,save_info.save_dir,fullfile(save_info.save_dir,'GCM',filesep));
    [gcm_dir,name,ext] = fileparts(save_info.nii_mat_name);
    try
        mkdir(gcm_dir)
    end
    save(save_info.nii_mat_name,'GCM');
end
if nargout==1
    varargout{1} = GCM;
end

