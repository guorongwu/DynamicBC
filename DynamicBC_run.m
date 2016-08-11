function DynamicBC_run(F)
%% batch of running DynamicBC
E = F.E;
R = cell2mat(get(E.rd_tvmodel,'val')); %sld,fls 
R1 = cell2mat(get(E.rd_rvw,'val')); %voxels ROI FCD 
R2 = cell2mat(get(E.rd_mr,'val')); %Default mask, User-Defined Mask,Nifti Label,TXT,Mat
flag_seed_ROI = cell2mat(get(E.rd_seed_ROI,'val'));

DataProcessDir = strcat(get(E.ed_bg_data_dir,'string'));
save_dir = get(E.ed_bg_dir_out,'string');
save_info.save_dir = save_dir;
prefix_name = get(E.ed_prefix(1),'string'); 
prefix_name_flag = length(strcat(prefix_name));
% log_name = get(E.ed_log_name(1),'string');
fc_ec = E.fc_ec ;

if fc_ec(1)
    str_fc_ec = 'FCM';
else
    str_fc_ec = 'GCM';
end

try
    mkdir(save_dir);
end

% global TV model parameter
if R(1) % sliding windows
    slw_id = find(cell2mat(get(E.rd_slw,'val'))>0);
    slw_alignment = get(E.rd_slw(slw_id),'string');
    slw_winsize = str2num(get(E.ed_rd11_tvmodel(1),'string'));
    slw_overlapsize = str2num(get(E.ed_rd11_tvmodel(2),'string'));
    slw_order = str2num(get(E.ed_rd11_tvmodel(3),'string'));
    save_info.slw_alignment = slw_id;    
    pvalue = str2num(get(E.ed_slw_pvalue,'string'));
    fwe_flag = get(E.rd_slw_conn(1),'val');
    
else %fls case
    fls_snr_id = find(cell2mat(get(E.rd_fls_snr,'val'))>0);
    mu = str2num(get(E.tx_fls_snr(fls_snr_id),'string'));
end

%% FC information
if any([R1([1 3]); (R1(2)&R2(3))])% seed FC/ FCD /ROIwise-Label
    disp('Default mode: fMRI(reading from NIFTI image), if not? choose "Set ROI/ ROI wise" ')
    sub = dir(DataProcessDir);
    sub(1:2)=[];
    k=[]; % in case non-directory file
    for i=1:length(sub)
        if ~sub(i).isdir
            k = [k i];
        end
    end
    sub(k)=[];    
    SubjectNum = length(sub);
    for i=1:SubjectNum
        SubjectID{i}=sub(i).name;
    end
else % ROIwise, i.e. mat/text file, for further auto select the mat/text file, if only select the directory
% %     sub = dir(DataProcessDir);
% %     SubjectNum = length(sub);
% %     for i=1:SubjectNum
% %         SubjectID{i}=sub(i).name;
% %     end
end

save_info.flag_par = matlabpool('size');

if any([R1([1 3]); (R1(2)&R2(3))]) % seed FC/ FCD
    %% read mask/label file. 
    if R2(1) % default mask
        fprintf('Generate Human brain mask from SPM default mask...\n')
        default_mask = fullfile(spm('Dir'),'apriori','brainmask.nii');
        if ~exist(default_mask)
            default_mask = fullfile(spm('Dir'),'tpm','brainmask.nii');
            if ~exist(default_mask)
                error('No human brain mask could be found! SPM should be installed.')
            end
        end        
        sub_dir = strcat(fullfile(DataProcessDir,SubjectID{1}));
        niiname = dir(fullfile(sub_dir,'*.nii'));
        if isempty(niiname)
            niiname = dir(fullfile(sub_dir,'*.img'));
            if isempty(niiname)             
                error(['no NIFTI file exist in: ',sub_dir]);
            end
        end
        niifile0 = fullfile(sub_dir,niiname(1).name);
        v0 = spm_vol(niifile0);
        mask_label_template = fullfile(save_dir,'brain_mask.nii');
        wgr_Reslice(default_mask,mask_label_template,v0(1).dim,0,niifile0);  
        fprintf('Brain mask: %s\n', mask_label_template)
    else
        mask_label_template = get(E.ed_mrs_files,'string');
    end
    
    vm = spm_vol(mask_label_template);
    mask_label_dat = spm_read_vols(vm);
    save_info.v = vm ; %% nifti header information
    save_info.v.dt = [16,0];
    save_info.v.n = [1 1];
    mask_label_dat(isnan(mask_label_dat))=0;
    disp('Default value 0/NaN is not in the mask/label!')

    if R1(1)
        if flag_seed_ROI(1)
            seed_ROI_mask = get(E.rd_seed_ROI(1),'TooltipString');
            seed_ROI_data = spm_read_vols(spm_vol(seed_ROI_mask));
            seed_ROI_data(isnan(seed_ROI_data))=0;
            save_info.seed_index = find(seed_ROI_data>0);
        else
            mni_radius = get(E.ed_MNI,'string') ;
            mni = [str2num(mni_radius{1}) str2num(mni_radius{2}) str2num(mni_radius{3})]; 
            radius = str2num(mni_radius{4});
            cor = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1)]*(inv(vm.mat))' ;
            cor(:,4) = [];
            cor = round(cor);
            [seed_cor] = wgr_center2ball(vm.dim,abs(diag(vm.mat(1:3,1:3)))',cor, radius) ; %ball
            seed_ind = sub2ind(vm.dim, seed_cor(:,1), seed_cor(:,2), seed_cor(:,3)) ;
            save_info.seed_index = seed_ind ;
        end
    else
        save_info.seed_index = [] ;
    end

    if any(R2([1 2])) % defaut mask/ specific mask.
        index = cell(1,1);
        index{1} = find(mask_label_dat>0);
        save_info.index = index{1};  %record mask index.
        nvar = length(index{1});
        disp('Default value>0 is inside the mask!');
        fprintf('There are %5.0f voxels inside the mask\n',nvar);

    elseif R2(3) % nifti label
        
        label_dat_id = unique(mask_label_dat(:)); 
        label_dat_id(label_dat_id==0)=[];
        nvar = length(label_dat_id);
        index = cell(nvar,1);
        for j = 1:nvar
            index{j}=find(mask_label_dat==label_dat_id(j));         
        end
        save_info.index = label_dat_id ; %record label id.
    end
    
    if R(1) % sliding windows
        if fwe_flag
            if fc_ec(1)
                pvalue = pvalue./(nvar*(nvar-1)/2); %FC
            else
                pvalue = pvalue./(nvar*(nvar-1));
            end
        end
    end
    
    %% read subject's data, for loop
    for isub=1:SubjectNum
        fprintf('Running subject %4.0f (all %4.0f subjects)\n',isub,SubjectNum);
        sub_dir = strcat(fullfile(DataProcessDir,SubjectID{isub}));
        niiname = dir(fullfile(sub_dir,'*.nii'));
        if isempty(niiname)
            niiname = dir(fullfile(sub_dir,'*.img'));
            if isempty(niiname)             
                warning(['Warning: no NIFTI file exist in: ',sub_dir]);
                continue;
            end
        end

        num_volume = length(niiname) ;
        if num_volume==1
            disp(['Only one file in: ',sub_dir])
            niifile = fullfile(sub_dir,niiname.name);
            v = spm_vol(niifile);
            data = spm_read_vols(v);
            nobs = size(data,4);
            if nobs<10
                disp(['Only ' num2str(nobs), ' time point in : ',niifile])
                disp(['We do not processed it, skip to next file.'])
                continue;
            else
                data = reshape(data,[],nobs); %nvar*nobs
            end
        elseif num_volume<10

            nobs = num_volume;
            disp(['Only ' num2str(nobs), ' Volume in : ',sub_dir])
            disp(['We do not processed it, skip to next file.'])
            continue;
        else
            nobs = num_volume;
            tmp_dat = spm_read_vols(spm_vol(fullfile(sub_dir,niiname(1).name)));
            data = zeros(prod(size(tmp_dat)),nobs); %nvar*nobs
            for j=1:nobs
                niifile = fullfile(sub_dir,niiname(j).name);
                v = spm_vol(niifile);
                tmp_dat = spm_read_vols(v);
                data(:,j) = tmp_dat(:);
            end
        end
        
        
        if any(R2([1 2])) % defaut mask/ specific mask.
            ROI_sig = data(index{1},:)';
        else
            ROI_sig = zeros(nobs,nvar);
            for j=1:nvar
                ROI_sig(:,j)=dynamicBC_nanmean(data(index{j},:),1);
            end
        end
        if isempty(save_info.seed_index) %seed ROI signal
            save_info.seed_signal = [];
        else
            save_info.seed_signal = dynamicBC_nanmean(data(save_info.seed_index,:),1)';
        end
            
        %% judge voxelswise/FCD/ROIwise-label
        
        if R1(1) % voxels

            save_info.flag_1n = 1;
            save_info.flag_nii = 1 ;
%             try
%                 mkdir( fullfile(save_dir, SubjectID{isub}) )
%             end
            if prefix_name_flag
                save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub},[prefix_name,'_',SubjectID{isub},'_',str_fc_ec,'.nii']);
            else
                save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub},[SubjectID{isub},'_',str_fc_ec,'.nii']);
            end

        elseif R1(3) % FCD

            save_info.flag_1n = 0;
            save_info.seed_index = [] ;
            save_info.flag_nii = 1 ;
%             try
%                 mkdir( fullfile(save_dir, SubjectID{isub}) )
%             end
            if prefix_name_flag
                save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub},[prefix_name,'_',SubjectID{isub},'_',str_fc_ec,'.nii']);
            else
                save_info.nii_mat_name = fullfile(save_dir,SubjectID{isub}, [SubjectID{isub},'_',str_fc_ec,'.nii']);
            end
            
        elseif R1(2)&R2(3) % ROIwise-label

            save_info.flag_1n = 0;
            save_info.seed_index = [] ;
            save_info.flag_nii = 0 ;
%             try
%                 mkdir( fullfile(save_dir, SubjectID{isub}) )
%             end
            if prefix_name_flag
                save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub},[prefix_name,'_',SubjectID{isub},'_',str_fc_ec,'.mat']);
            else
                save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub},[SubjectID{isub},'_',str_fc_ec,'.mat']);
            end

        end
        
        if R(1) % sliding windows
            if fc_ec(1)
                DynamicBC_sliding_window_FC(ROI_sig,slw_winsize,slw_overlapsize,pvalue,save_info);
            elseif fc_ec(2)
                DynamicBC_sliding_window_GC(ROI_sig,slw_winsize,slw_overlapsize,pvalue,slw_order,save_info);
            end
        else % FLS
            DynamicBC_fls_FC(ROI_sig,mu,save_info);
        end
    end
    
else % ROI-wise, i.e. txt or mat file
    
    tx_mat_file = get(E.ed_mrs_files,'string');
    save_info.flag_1n = 0;
    save_info.seed_index = [] ;
    save_info.flag_nii = 0 ;
    SubjectNum = size(tx_mat_file,1);
    SubID = cell(SubjectNum,2);
    for isub=1:SubjectNum
        fprintf('Running subject %4.0f (all %4.0f subjects)\n',isub,SubjectNum);
        subj_file = strcat(tx_mat_file(isub,:));
        [sub_path, SubjectID{isub}, ext]= fileparts(subj_file);
        if isub>1
           TF = strcmpi(SubjectID{isub},SubjectID(1:isub-1));
           if sum(TF(:))>0
               SubjectID{isub} = [SubjectID{isub},'_A'];
           end
        end
        SubID{isub,1} = subj_file;
        SubID{isub,2} = SubjectID{isub};
        
        if R2(4)
            ROI_sig = load(subj_file);
        elseif R2(5)
            load(subj_file,get(E.ed_mat,'string'));
            eval(['ROI_sig = ',get(E.ed_mat,'string'),';']);
        end
        nvar = size(ROI_sig,2);
        
        if R(1) % sliding windows
            if fwe_flag
                if fc_ec(1)
                    pvalue = pvalue./(nvar*(nvar-1)/2); %FC
                else
                    pvalue = pvalue./(nvar*(nvar-1));
                end
            end
        end
        
        if nvar>1500
            error('FLS do not support number of variabel exceed 1500!(cost too much time & memory)')
        end
        
%         try
%             mkdir( fullfile(save_dir, SubjectID{isub}) )
%         end
        if prefix_name_flag
            save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub}, [prefix_name,'_',SubjectID{isub},'_',str_fc_ec,'.mat']);
        else
            save_info.nii_mat_name = fullfile(save_dir, SubjectID{isub}, [SubjectID{isub},'_',str_fc_ec,'.mat']);
        end

        if R(1) % sliding windows
            if fc_ec(1)
                DynamicBC_sliding_window_FC(ROI_sig,slw_winsize,slw_overlapsize,pvalue,save_info);
            elseif fc_ec(2)
                DynamicBC_sliding_window_GC(ROI_sig,slw_winsize,slw_overlapsize,pvalue,slw_order,save_info);
            end                
        else % FLS
            DynamicBC_fls_FC(ROI_sig,mu,save_info);
        end
    end
    xlswrite(fullfile(save_dir,'Subject_ID_information_txt_mat.xls'),SubID);
end

set(F.pb_run_check(1),'backgroundc',get(F.hfig,'color'));
fprintf('====Finish DynamicBC!^_^====\n\n')
delete(F.hfig);

%% subfunction

function wgr_Reslice(PI,PO,NewVoxSize,hld,TargetSpace)
%   PI - input filename
%   PO - output filename
%   NewVoxSize - 1x3 matrix of new vox size.
%   hld - interpolation method. 0: Nearest Neighbour. 1: Trilinear.
%   TargetSpace - Define the target space. 'ImageItself': defined by the image itself (corresponds  to the new voxel size); 'XXX.img': defined by a target image 'XXX.img' (the NewVoxSize parameter will be discarded in such a case).
if nargin<=4
    TargetSpace='ImageItself';
end

if ~strcmpi(TargetSpace,'ImageItself')
    headIN = spm_vol(TargetSpace) ;
    headIN = headIN(1);
    dataIN = spm_read_vols(headIN);
    mat=headIN.mat;
    dim=headIN.dim;
else
    headIN = spm_vol(PI) ;
    headIN = headIN(1);
    dataIN = spm_read_vols(headIN);
    origin=headIN.mat(1:3,4);
    origin=origin+[headIN.mat(1,1);headIN.mat(2,2);headIN.mat(3,3)]-[NewVoxSize(1)*sign(headIN.mat(1,1));NewVoxSize(2)*sign(headIN.mat(2,2));NewVoxSize(3)*sign(headIN.mat(3,3))];
    origin=round(origin./NewVoxSize').*NewVoxSize';
    mat = [NewVoxSize(1)*sign(headIN.mat(1,1))                 0                                   0                       origin(1)
        0                         NewVoxSize(2)*sign(headIN.mat(2,2))              0                       origin(2)
        0                                      0                      NewVoxSize(3)*sign(headIN.mat(3,3))  origin(3)
        0                                      0                                   0                          1      ];

    dim=(headIN.dim-1).*diag(headIN.mat(1:3,1:3))';
    dim=floor(abs(dim./NewVoxSize))+1;
end
VI          = spm_vol(PI);
VI   = VI(1);
VO          = VI;
VO.fname    = deblank(PO);
VO.mat      = mat;
VO.dim(1:3) = dim;

[x1,x2] = ndgrid(1:dim(1),1:dim(2));
d     = [hld*[1 1 1]' [1 1 0]'];
C = spm_bsplinc(VI, d);
v = zeros(dim);
for x3 = 1:dim(3),
    [tmp,y1,y2,y3] = getmask(inv(mat\VI.mat),x1,x2,x3,VI.dim(1:3),[1 1 0]');
    v(:,:,x3)      = spm_bsplins(C, y1,y2,y3, d);
end;
VO = spm_write_vol(VO,v);
%_______________________________________________________________________
function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = logical(ones(size(y1)));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end;
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end;
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end;
return;
%_______________________________________________________________________

function [cor_s] = wgr_center2ball(ABrainSize, AVoxelSize, AROICenter, AROIRadius);
%ABrainSize, such as [61, 73, 61]
%AROICenter ,cor coordinate
radiusX =round(AROIRadius /AVoxelSize(1));
if (AROICenter(1)-radiusX)>=1 && (AROICenter(1)+radiusX)<=ABrainSize(1)
    rangeX	=(AROICenter(1)-radiusX):(AROICenter(1)+radiusX);
elseif (AROICenter(1)-radiusX)<1 && (AROICenter(1)+radiusX)<=ABrainSize(1)
    rangeX	=1:(AROICenter(1)+radiusX);
elseif (AROICenter(1)-radiusX)>=1 && (AROICenter(1)+radiusX)>ABrainSize(1)
    rangeX	=(AROICenter(1)-radiusX):ABrainSize(1);
else
    rangeX =1:ABrainSize(1);
end

radiusY =round(AROIRadius /AVoxelSize(2));
if (AROICenter(2)-radiusY)>=1 && (AROICenter(2)+radiusY)<=ABrainSize(2)
    rangeY	=(AROICenter(2)-radiusY):(AROICenter(2)+radiusY);
elseif (AROICenter(2)-radiusY)<1 && (AROICenter(2)+radiusY)<=ABrainSize(2)
    rangeY	=1:(AROICenter(2)+radiusY);
elseif (AROICenter(2)-radiusY)>=1 && (AROICenter(2)+radiusY)>ABrainSize(2)
    rangeY	=(AROICenter(2)-radiusY):ABrainSize(2);
else
    rangeY =1:ABrainSize(2);
end

radiusZ =round(AROIRadius /AVoxelSize(3));
if (AROICenter(3)-radiusZ)>=1 && (AROICenter(3)+radiusZ)<=ABrainSize(3)
    rangeZ	=(AROICenter(3)-radiusZ):(AROICenter(3)+radiusZ);
elseif (AROICenter(3)-radiusZ)<1 && (AROICenter(3)+radiusZ)<=ABrainSize(3)
    rangeZ	=1:(AROICenter(3)+radiusZ);
elseif (AROICenter(3)-radiusZ)>=1 && (AROICenter(3)+radiusZ)>ABrainSize(3)
    rangeZ	=(AROICenter(3)-radiusZ):ABrainSize(3);
else
    rangeZ =1:ABrainSize(3);
end
mni_s = []; cor_s = [];
for x=rangeX, for y=rangeY, for z=rangeZ,
    %Ball Definition, Computing within a cubic to minimize the time to be consumed
    if norm(([x, y, z] -AROICenter).*AVoxelSize)<=AROIRadius,
        cor_s = [cor_s; [x, y, z]];
    end
end; end; end;

function m = dynamicBC_nanmean(x,dim)
%NANMEAN Mean value, ignoring NaNs.
%   M = NANMEAN(X) returns the sample mean of X, treating NaNs as missing
%   values.  For vector input, M is the mean value of the non-NaN elements
%   in X.  For matrix input, M is a row vector containing the mean value of
%   non-NaN elements in each column.  For N-D arrays, NANMEAN operates
%   along the first non-singleton dimension.
%
%   NANMEAN(X,DIM) takes the mean along dimension DIM of X.
%
%   See also MEAN, NANMEDIAN, NANSTD, NANVAR, NANMIN, NANMAX, NANSUM.

%   Copyright 1993-2004 The MathWorks, Inc.


% Find NaNs and set them to zero
nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
end
