function [varargout] = DynamicBC_sliding_window_FC(data,window,overlap,pvalue,save_info)
%% calculate sliding window functional connectivity (bivariate)
% overlap = 0.1; %e.g. time bin: [1:50],[46:95], [91:140],...
% window = 50;
[nobs, nvar] = size(data);
step=ceil(window-overlap*window); % 10% overlap
% step=ceil((1-overlap)*window); % 10% overlap
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
    slides=floor((nobs-window)/step)+1;
end
num0 = ceil(log10(slides))+2;
FC = cell(slides,1);  
t1=1-step;
t2=window-step;
%sliding window 
nii_name = cell(slides,4);

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
    if save_info.flag_nii % save nii: seed FC,FCD.
        v = save_info.v;
        v.fname = strcat(save_info.nii_mat_name);
        [pathstr, name, ext] = fileparts(v.fname);       
       
        num1 = num0 - length(num2str(k1));
        data_save = zeros(save_info.v.dim);
        if ~save_info.flag_1n % FCD
            pathstr_fcd = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FCD_SW_map',filesep));
            if k==1
                try
                    mkdir(pathstr_fcd)
                end
            end
            [FCM] = wgr_correl_FCD(dat,nvar, pvalue);
            v.fname = fullfile(pathstr_fcd,[name,'_abs_wei_',repmat('0',1,num1),num2str(k1),ext]);
            data_save(save_info.index) = FCM.fcd_abs_wei; 
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;
            
            v.fname = fullfile(pathstr_fcd,[name,'_abs_bin_',repmat('0',1,num1),num2str(k1),ext]);
            data_save(save_info.index) = FCM.fcd_abs_bin ;
            spm_write_vol(v,data_save);
            nii_name{k,2} = v;
            
            v.fname = fullfile(pathstr_fcd,[name,'_pos_wei_',repmat('0',1,num1),num2str(k1),ext]);
            data_save(save_info.index) = FCM.fcd_pos_wei ;
            spm_write_vol(v,data_save);
            nii_name{k,3} = v;
            
            v.fname = fullfile(pathstr_fcd,[name,'_pos_bin_',repmat('0',1,num1),num2str(k1),ext]);
            data_save(save_info.index) = FCM.fcd_pos_bin ;
            spm_write_vol(v,data_save);
            nii_name{k,4} = v;
            
        else  %if save_info.flag_1n % seed FC
            
            FCM = wgr_correl(dat,nvar,save_info.flag_1n,pvalue,save_info.seed_signal(t1:t2,:));
            Matrix = zeros(FCM.nvar,1);
            for i=1:FCM.nbin
                Matrix(FCM.indX{i},1) =  FCM.matrix{i,1};
            end
            data_save(save_info.index) = Matrix;
            pathstr_r = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'seed_CORR_FCmap',filesep));
            if k==1
                try
                    mkdir(pathstr_r)
                end
            end
            v.fname = fullfile(pathstr_r,[name,'_',repmat('0',1,num1),num2str(k1),ext]);
            spm_write_vol(v,data_save);
            nii_name{k,1} = v;
            
            pathstr_z = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'seed_Z_FCmap',filesep));
            if k==1
                try
                    mkdir(pathstr_z)
                end
            end
            v.fname = fullfile(pathstr_z,['Z_',name,'_',repmat('0',1,num1),num2str(k1),ext]);
            data_save(save_info.index) = atanh(Matrix); %0.5 * log((1 + Matrix)./(1 - Matrix));
            spm_write_vol(v,data_save);
            nii_name{k,2} = v;
        end
        

    else % ROIwise save as mat file
        
        FCM = wgr_correl(dat,nvar,save_info.flag_1n,pvalue,[]);
        Matrix = sparse(FCM.nvar,FCM.nvar);
        for i=1:FCM.nbin
            for j=1:FCM.nbin
                Matrix(FCM.indX{i},FCM.indX{j}) =  FCM.matrix{i,j};
            end
        end
        FC{k,1}  = Matrix; 
        time_alignment(k) = k1; 
    end
end

%% variance calculation

if save_info.flag_nii % save nii: seed FC,FCD.
    if ~save_info.flag_1n % FCD
        pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FCD_SW_Variance',filesep));
        try
            mkdir(pathstr_v)
        end
        str_typ = {'_abs_wei_variance','_abs_bin_variance','_pos_wei_variance','_pos_bin_variance'};
        str_typ2 = {'_abs_wei_std','_abs_bin_std','_pos_wei_std','_pos_bin_std'};
        str_typ3 = {'_abs_wei_mean','_abs_bin_mean','_pos_wei_mean','_pos_bin_mean'};
        str_typ4 = {'_abs_wei_CV','_abs_bin_CV','_pos_wei_CV','_pos_bin_CV'};
        for j=1:4
            data0 = zeros(slides, nvar);
            for k=1:slides
                v = spm_vol(nii_name{k,j}.fname);
                dat = spm_read_vols(v);
                data0(k,:) = dat(save_info.index);
            end
            datav = var(data0,0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ{j},ext]);
            data_save(save_info.index) = datav; 
            spm_write_vol(v,data_save);
            
            datas = std(data0,0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ2{j},ext]);
            data_save(save_info.index) = datas; 
            spm_write_vol(v,data_save);
            
            datam = mean(data0,1);
            v.fname = fullfile(pathstr_v,[name,str_typ3{j},ext]);
            data_save(save_info.index) = datam; 
            spm_write_vol(v,data_save);
            
            datacv = datas./datam;
            v.fname = fullfile(pathstr_v,[name,str_typ4{j},ext]);
            data_save(save_info.index) = datacv; 
            spm_write_vol(v,data_save);

            v.fname = fullfile(pathstr_v,[name,str_typ4{j},'_abs',ext]);
            data_save(save_info.index) = abs(datacv); 
            spm_write_vol(v,data_save);
            
        end
    else
        if slides>1 
            pathstr_v = strrep(pathstr,save_info.save_dir,fullfile(save_info.save_dir,'FC_SW_Variance',filesep));
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
                    prefix = '';
                else
                    prefix =  'Z_';
                end
                data = var(data0,0,1);
                v.fname = fullfile(pathstr_v,[prefix,name,'_variance',ext]);
                data_save(save_info.index) = data; 
                spm_write_vol(v,data_save);
                datas = std(data0,0,1);
                v.fname = fullfile(pathstr_v,[prefix,name,'_std',ext]);
                data_save(save_info.index) = datas; 
                spm_write_vol(v,data_save);
                datam = mean(data0,1);
                v.fname = fullfile(pathstr_v,[prefix,name,'_mean',ext]);
                data_save(save_info.index) = datam; 
                spm_write_vol(v,data_save);
                cv = datas./datam;
                v.fname = fullfile(pathstr_v,[prefix,name,'_CV',ext]);
                data_save(save_info.index) = cv; 
                spm_write_vol(v,data_save);
                v.fname = fullfile(pathstr_v,[prefix,name,'_CV_abs',ext]);
                data_save(save_info.index) = abs(cv); 
                spm_write_vol(v,data_save);
            end
        end
    end
    
else
    
    data_var = zeros(slides, nvar*nvar);
    for k=1:slides
        tmp_matrix  = FC{k,1};
        data_var(k,:) = tmp_matrix(:); 
    end
    data_var = var(data_var,0,1);
    data_var = reshape(data_var,nvar,nvar);
    
end

if ~save_info.flag_nii
    r_th = FCM.r_th;
    clear FCM
    FCM.pvalue = pvalue;
    FCM.r_th = r_th;
    FCM.Matrix = FC;
    FCM.time_alignment = time_alignment;
    FCM.variance = data_var;
    save_info.nii_mat_name = strrep(save_info.nii_mat_name,save_info.save_dir,fullfile(save_info.save_dir,'FCM',filesep));
    [fcm_dir] = fileparts(save_info.nii_mat_name);
    try
        mkdir(fcm_dir)
    end
    save(save_info.nii_mat_name,'FCM');
end
if nargout==1
    varargout{1} = FCM;
end
clear data0

function [M] = wgr_correl(data,nvar,flag_1n, pvalue,seed_signal,m)
% Compute correlation matrix A.
if nargin<6
    m = 5000; %block size
end
nobs = size(data,1);
M.nvar=nvar;
M.nobs=nobs;
M.pvalue=pvalue;
r_th = invpearson(pvalue,nobs);
M.r_th=r_th;

if nvar <= m 
    M.nbin = 1;
    M.indX{1} = 1:nvar;
    if flag_1n %seed FC  1 x n
        A = corr(seed_signal,data);
        M.matrix{1,1} = A ;
    else % n x n
        A = corr(data);
        M.matrix{1,1} = A ;
    end       
    
else 
    
    nbin = ceil(nvar/m);
    M.nbin = nbin;
    for i=1:nbin
        if i~=nbin
            ind_X = (i-1)*m+1:i*m ;   
        else
            ind_X = (i-1)*m+1:nvar ;
        end

        X=data(:,ind_X);
        M.indX{i} = ind_X; 
        
        if flag_1n %seed FC: 1 x n
            RR = corr(seed_signal,X);
            M.matrix{i,1} = RR;
        else  % FC: n x n
            for j=1:nbin
                if j~=nbin
                    ind_Y = (j-1)*m+1:j*m ;
                else
                    ind_Y = (j-1)*m+1:nvar ;
                end

                Y=data(:,ind_Y);

                RR = corr(X, Y);        

                M.matrix{i,j} = RR;
%                 if flag_wei
%                     M.matrix{i,j} = RR.*(abs(RR)>r_th);
%                 else
%                     M.matrix{i,j} = RR ;
%                 end
            end
        end
    end  
end
    
function [M] = wgr_correl_FCD(data,nvar, pvalue)
% Compute correlation matrix A.
m = 5000; %block size
nobs = size(data,1);
M.nvar=nvar;
M.nobs=nobs;

M.pvalue=pvalue;
r_th = invpearson(pvalue,nobs);
M.r_th=r_th;

if nvar <= m 
    M.nbin = 1;
    M.indX{1} = 1:nvar;
    A = corr(data);
    A(isnan(A)) = 0;
    comm = {'abs(A); MA = MA.*(MA>r_th); MA = atanh(MA); MA(isinf(MA)) = 0;','A.*(A>r_th); MA = atanh(MA); MA(isinf(MA)) = 0;','abs(A)>r_th','A>r_th'};
    comm_str = {'fcd_abs_wei','fcd_pos_wei','fcd_abs_bin','fcd_pos_bin'};
    MA = zeros(nvar);
    for i=1:length(comm)
        if i<3
            str_commond = ['MA = ',comm{i},'; M.',comm_str{i},' = sum(MA,2);'];
        else
            str_commond = ['MA = ',comm{i},'; M.',comm_str{i},' = sum(MA,2)-1;'];
        end
%         disp(str_commond)
        eval(str_commond);
    end
else 
    nbin = ceil(nvar/m);
    M.nbin = nbin;
    for i=1:nbin
        if i~=nbin
            ind_X = (i-1)*m+1:i*m ;   
        else
            ind_X = (i-1)*m+1:nvar ;
        end

        X=data(:,ind_X);
        M.indX{i} = ind_X; 
        matrixs = zeros(length(ind_X),nvar);
        for j=1:nbin
            if j~=nbin
                ind_Y = (j-1)*m+1:j*m ;
            else
                ind_Y = (j-1)*m+1:nvar ;
            end
            
            Y=data(:,ind_Y);            
            RR = corr(X, Y);
            RR(isnan(RR)) = 0;
            matrixs(:,ind_Y) = RR;
        end
        
        comm = {'abs(matrixs);MA = MA.*(MA>r_th); MA = atanh(MA); MA(isinf(MA)) = 0','matrixs.*(matrixs>r_th); MA = atanh(MA); MA(isinf(MA)) = 0','abs(matrixs)>r_th','matrixs>r_th'};
        comm_str = {'fcd_abs_wei','fcd_pos_wei','fcd_abs_bin','fcd_pos_bin'};
        MA = zeros(length(ind_X),nvar);
        for i0=1:length(comm)
            if i0<3
                str_commond = ['MA = ',comm{i0},'; M.',comm_str{i0},'(ind_X,1) = sum(MA,2);'];
            else
                str_commond = ['MA = ',comm{i0},'; M.',comm_str{i0},'(ind_X,1) = sum(MA,2)-1;'];
            end
%             disp(str_commond)
            eval(str_commond);
        end
       
    end
end

function [r] = invpearson(p,n)
% from p-value to Pearson's linear correlation value . 
t = tinv ( 1-p/2, n-2) ;
r = sqrt(t.^2./(n-2 + t.^2)) ;    