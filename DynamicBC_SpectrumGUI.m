function DynamicBC_SpectrumGUI()
%version 2.0 2018.10.20
D.fig = figure('Name','Amplitude Spectrum (v2.0)',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'unit','normalized',...
    'color',[0.95 0.95 0.95],...
    'position',[0.25 0.2 0.6 0.25]);
movegui(D.fig,'center');           
            
%%=========================Data information=====================================
pos_Data = [0.025 0.025 0.7 0.95];
pos_vis = [0.75 0.025 0.23 0.95];

D.data_info = uibuttongroup('units','norm',...
    'Title','Data Information','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
                     'pos',pos_Data);   

pos_data_type = [0.05 0.71 0.39 0.2];
pos_FS_str = [0.45 0.71 0.35 0.2];
pos_FS_input = [0.8 0.71 0.15 0.2];
pos_data_mask_str  =  [0.05 0.49 0.39 0.2];
pos_data_mask_input  =  [0.45 0.49 0.5 0.2];
pos_var_str = [0.05 0.49 0.25 0.2];
pos_var_input =  [0.31 0.49 0.2 0.2];
pos_connind_str = [0.52 0.5 0.28 0.2] ;
pos_connind_input = [0.8 0.5 0.15 0.2] ;
pos_data_select_str = [0.05 0.27 0.35 0.2];
pos_data_select = [0.4 0.27 0.55 0.2];
pos_out_str = [0.05 0.05 0.35 0.2];
pos_out_select = [0.4 0.05 0.55 0.2];

                
D.data_type = uicontrol('Parent',D.data_info,...
                    'style','popup',...
                    'units','norm',...
                    'position',pos_data_type,...
                    'string',{'NIFTI image','Matrix','Image or Matrix?'},...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...                    
                    'horizontalalign','center',...                    
                    'value',3,...
                    'TooltipString','...');  
D.tx_fs = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_FS_str,...
                    'string','Sampling Rate:',...                                      
                    'fontunits', 'normalized', 'enable','off','fontsize',0.5,...                   
                    'horizontalalign','Right');  
D.ed_fs = uicontrol('Parent',D.data_info,...
                    'style','edit',...
                    'units','norm',...
                    'position',pos_FS_input,...
                    'string','1/2',...                   
                    'fontunits', 'normalized', 'fontsize',0.5,...              
                    'horizontalalign','center',...                   
                    'visible','on',...
                    'TooltipString','1/TR');   
% NIFTI image
D.tx_mask = uicontrol('Parent',D.data_info,...
                'style','popup',...
                'units','norm',...
                'position',pos_data_mask_str,...
                'string',{'Mask(NIFTI)','Sphere([x,y,z,radius])','Mask(NIFTI/MNI Coordinates)?'},...
                'val',3,...
                'backgroundc',get(D.data_info,'backgroundc'),...                     
                'visible','off',...
                'TooltipString','Click me to select NIFTI image',...
                'fontunits', 'normalized', 'fontsize',0.5,...            
                'horizontalalign','left');    
                
D.ed_mask = uicontrol('Parent',D.data_info,...
                        'style','edit',...
                        'unit','norm',...
                        'position',pos_data_mask_input,...
                        'fontunits', 'normalized', 'fontsize',0.5,...                         
                        'visible','off',...
                        'string','',...
                        'horizontalalign','right');  
              
%% Matrix                        
D.tx_varname = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_var_str,...
                    'string','Variable:',...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...                   
                    'value',0,...                     
                    'visible','off',...                   
                    'TooltipString','e.g. FCM.Matrix');                 
D.ed_varname = uicontrol('Parent',D.data_info,...
                    'style','edit',...
                    'units','norm',...
                    'position',pos_var_input,...
                    'string','FCM.Matrix',...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...               
                    'horizontalalign','center',...                   
                    'visible','off',...
                    'TooltipString','FCM.Matrix or GCM.Matrix');                            
   
 D.tx_ID = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_connind_str,...
                    'string','Matrix index:',...                    
                    'fontunits', 'normalized', 'fontsize',0.5,...                 
                    'visible','off',...                       
                    'value',0);                 
D.ed_ID= uicontrol('Parent',D.data_info,...
                    'style','edit',...
                    'units','norm',...
                    'position',pos_connind_input,...
                    'string','[1, 2]',...                    
                    'fontunits', 'normalized', 'fontsize',0.5,...                   
                    'horizontalalign','center',...                   
                    'visible','off',...
                    'TooltipString','e.g. [2,3]: row 2,column 3');     
                       
%% Data selection
D.data_tx = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_data_select_str,...
                    'string','Data Folder:',...                
                    'fontunits', 'normalized', 'fontsize',0.5,...                  
                    'horizontalalign','center',...                    
                    'value',1,...
                    'TooltipString','...');    
D.ed_data = uicontrol(D.data_info,...
                    'style','edit',...
                    'unit','norm',...
                    'position',pos_data_select,...
                    'fontunits', 'normalized', 'fontsize',0.5,...
                    'string','not selected');    
              
%% OutputDirectory
D.out_str = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_out_str,...
                    'string','OUT Folder:',...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...                   
                    'horizontalalign','center',...                   
                    'value',1,...
                    'TooltipString','...');    
D.ed_outdir = uicontrol(D.data_info,...
                    'style','edit',...
                    'unit','norm',...
                    'position',pos_out_select,...
                    'fontunits', 'normalized', 'fontsize',0.5,...
                    'string','not selected');              

% vis
pos_sub = [0.05 0.75 0.9 0.2] ;
pos_matrix = [0.05 0.5 0.9 0.2] ;
pos_plot = [0.05 0.05 0.9 0.2] ;
pos_Run = [0.05 0.27 0.9 0.2] ;

D.vis = uibuttongroup('Parent',D.fig,...
                     'units','norm',...
                     'Title','Subject/Matrix','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[1 1 1]*0.5,...
                     'pos',pos_vis);   
D.sub_select = uicontrol('Parent',D.vis,...
                'style','popup',...
                'units','norm',...
                'position',pos_sub,...
                'string',{'All subjects','Select parts'},...               
                'fontunits', 'normalized', 'fontsize',0.5,...                
                'foregroundcolor',[.1 .1 .1],...
                'TooltipString','Choose some subjects',...
                'value',1); 

D.matrix_select = uicontrol('Parent',D.vis,...
                'style','popup',...
                'units','norm',...
                'position',pos_matrix,...
                'string',{'Full Matrix','Only selected matrix index'},...               
                'fontunits', 'normalized', 'fontsize',0.5,...                
                'foregroundcolor',[.1 .1 .1],...
                'TooltipString','perform calculation on all/one elements of the matrix',...
                'vis','off',...
                'value',2);     

D.FCGC_select = uicontrol('Parent',D.vis,...
                'style','popup',...
                'units','norm',...
                'position',pos_matrix,...
                'string',{'*','IN','OUT'},...               
                'fontunits', 'normalized', 'fontsize',0.5,...                
                'foregroundcolor',[.1 .1 .1],...
                'TooltipString','perform calculation on all/one elements of the matrix',...
                'vis','off',...
                'value',1);  
            
D.pb_run = uicontrol(D.vis,...
                  'style','pushbutton',...
                  'units','norm',...
                  'position',pos_Run,...
                  'string','RUN',...                
                  'TooltipString','Click me!',...
                  'fontunits', 'normalized', 'fontsize',0.6,'val',1); 
D.pb_plot_fig = uicontrol(D.vis,...
                  'style','pushbutton',...
                  'units','norm',...
                  'position',pos_plot,...
                  'string','Results Plot',...
                  'TooltipString','Click me After RUN!',...
                  'fontunits', 'normalized', 'fontsize',0.6);   

set(D.tx_mask,'callback',{@wgr_select_dir_rewrite,D,'mask',''});
set(D.data_tx,'callback',{@wgr_select_dir_rewrite,D,'data',''});
set(D.out_str,'callback',{@wgr_select_dir_rewrite,D,'output',''});
set(D.data_type,'callback',{@wgr_data_type,D});
set(D.pb_run,'callback',{@wgr_run_call,D});
set(D.pb_plot_fig,'callback',{@wgr_plot_fig_call,D});
set(D.sub_select,'callback',{@wgr_sub_select_call,D});

function []=wgr_run_call(varargin)
D = varargin{3};  % Get structure.
Datafold = get(D.ed_data,'string');
SubFold = dir(Datafold); 
if isempty(SubFold)
    error(sprintf('No subject data in %s\n',Datafold))
    return;
else
    SubFold(1:2)=[];
    tmp=[];
    for i=1:length(SubFold) %check for mac os
        if ~SubFold(i).isdir
            tmp = [tmp i];
        end
    end
    SubFold(tmp)=[]; clear tmp
end
outputd = get(D.ed_outdir,'string');
Fs = str2num(get(D.ed_fs,'string'));
flag_datatype = get(D.data_type,'val');
flag_seed_ROI = get(D.tx_mask,'val');
varname = get(D.ed_varname,'string');
matfull = get(D.matrix_select,'val');
NII_ioflag = get(D.FCGC_select,'val');
iostr_filter0 = get(D.FCGC_select,'string');
iostr_filter = iostr_filter0{NII_ioflag};
id = str2num(get(D.ed_ID,'string'));

if flag_datatype==1
    sub01dir = fullfile(Datafold,SubFold(1).name);
    if NII_ioflag==1
        niifile = spm_select('FPList',sub01dir,'image');
    else
        niifile = spm_select('FPList',sub01dir,'image',[iostr_filter,'.*$']);
    end
    
    ROI_mask = get(D.ed_mask,'String');
    if flag_seed_ROI==2
        vm = spm_vol(strcat(niifile(1,:)));
        mni_radius = str2num(ROI_mask);
        mni = mni_radius(1:3);
        radius = mni_radius(4);
        %%
        cor = [mni(:,1) mni(:,2) mni(:,3) ones(size(mni,1),1) ]*(inv(vm.mat))' ;
        cor(:,4) = [];
        cor = round(cor);
        [seed_cor] = wgr_center2ball(vm.dim,abs(diag(vm.mat(1:3,1:3)))',cor, radius) ; %ball
        seed_ind = sub2ind(vm.dim, seed_cor(:,1), seed_cor(:,2), seed_cor(:,3)) ;
    elseif flag_seed_ROI==1
        
        niifile0 = strcat(niifile(1,:));
        v0 = spm_vol(niifile);
        vs = spm_vol(ROI_mask);
        if any(v0(1).dim-vs(1).dim)
            fprintf('%s, DIM:[%d,%d,%d]\n',niifile0,v0(1).dim)
            fprintf('Seed ROI mask %s, DIM:[%d,%d,%d]\n',ROI_mask,vs(1).dim)
            fprintf('(Nearest Neighbour)Reslice ROI mask ... \n')
            [fpath_seed,name,ext] = fileparts(ROI_mask);
            tmpext = strfind(ext,',');
            if ~isempty(tmpext)
                ext = ext(1:tmpext(1)-1);
            end
            ROI_mask0 = fullfile(fpath_seed,[name,'_resliced',ext]);
            DynamicBC_Reslice(ROI_mask,ROI_mask0,0,niifile0);
            vs = spm_vol(ROI_mask0);
        end
        ROI_data = spm_read_vols(vs);
        ROI_data(isnan(ROI_data))=0;
        seed_ind = find(ROI_data);
    end
end

subid = str2num(get(D.sub_select,'TooltipString')) ;
if isempty(subid)
    subid = 1:length(SubFold);
end

if ~exist(outputd,'dir')
    mkdir(outputd)
end
y={};  idnewf=1;
if flag_datatype==1
    for isub = subid
        subdir = fullfile(Datafold,SubFold(isub).name);
        if NII_ioflag==1
            niifile = spm_select('FPList',subdir,'image');
        else
            niifile = spm_select('FPList',subdir,'image',[iostr_filter,'.*$']);
        end
        data = spm_read_vols(spm_vol(niifile));
        data = reshape(data,[],size(data,4))';
        y{isub,1}= dynamicBC_nanmean(data(:,seed_ind),2);        
    end
    save(fullfile(outputd,'power_spec.mat'),'y','Fs','subid','Datafold','SubFold','idnewf','flag_datatype','ROI_mask','matfull','NII_ioflag')
else
    
    str = ['matrix = ',varname,';'];
    fprintf([str,'\n'])
    for isub = subid
        subdir = fullfile(Datafold,SubFold(isub).name);
        matfile = spm_select('FPList',subdir,'mat');
        if ~isempty(matfile)&&size(matfile,1)==1
            load(matfile)
            eval(str)
            nobs = length(matrix);
            if matfull==2
                y0 = zeros(nobs,1);
%                 fprintf('# of observation: %d\n',nobs)
                for j=1:nobs
                    tmp = matrix{j};
                    y0(j) = tmp(id(1),id(2));
                end
            else
                if isub==subid(1)
                    tt = matrix{1}-matrix{1}'; tt(tt<10^-6)=0;
                    info.dim = size(tt);
                    if any(tt(:))
                        idf = find([ones(info.dim)-eye(info.dim)]);
                        disp('asymmetry matrix')
                        info.symmetry = 0;
                        idnew = sub2ind(info.dim,id(1),id(2));
                    else
                        idf = find(tril(ones(size(tt)),-1));
                        disp('symmetry matrix')
                        info.symmetry = 1;
                        idnew = sub2ind(info.dim,max(id),min(id)); %lower triangular matrix
                    end
                    
                    idnewf = find(idf==idnew, 1);

                    if isempty(idnewf)
                        aa = sfprintf('matrix index error! data dimension is %dx%d, input ID is %d (=(%d-1)*%d+%d)\n',size(tt),idnew,id(2),size(tt,1),id(1));
                        error(aa);
                    end
                    clear tt  
                end
                y0=zeros(nobs,length(idf));
                for j=1:nobs
                    tmp = matrix{j};
                    y0(j,:) = tmp(idf);
                end
            end
        else
            y0=[];
        end
        y{isub,1}=y0;
    end
    save(fullfile(outputd,'power_spec.mat'),'y','Fs','subid','Datafold','SubFold','idnewf','idf','flag_datatype','matfull','info','id')
end
wgr_plot_powerspectral(y,Fs,flag_datatype,matfull,subid,SubFold,idnewf)

function wgr_plot_fig_call(varargin)
D = varargin{3};
Datafold0 = get(D.ed_data,'string');
outputd = get(D.ed_outdir,'string');
savemat = fullfile(outputd,'power_spec.mat');
fprintf('Load data: %s\n',savemat)
flag_datatype0 = get(D.data_type,'val');
ROI_mask0 = get(D.ed_mask,'string');
id0 = str2num(get(D.ed_ID,'string'));
matfull0 = get(D.matrix_select,'val');
NII_ioflag0 = get(D.FCGC_select,'val');
iostr_filterc = get(D.FCGC_select,'string');

subid0 = str2num(get(D.sub_select,'TooltipString')) ;
load(savemat)
if isempty(subid0)
    subid0 = 1:length(SubFold);
end
if flag_datatype==1
    if NII_ioflag0~=NII_ioflag
        aaa=sprintf('different Data filter! please check Subject/Matrix information pannel\n;%s\n %s\n',iostr_filterc{NII_ioflag0},iostr_filterc{NII_ioflag});
        error(aaa)
    end
end
if flag_datatype==2&&any(id0-id)
    if matfull==2
        fprintf('Click Run then plot \n')
    else
        if info.symmetry
            try
                idnew = sub2ind(info.dim,max(id0),min(id0)); %lower triangular matrix
            catch
                ff = sprintf('please reset matrix index, maximum is [%d %d] \n',info.dim);
                error(ff)
            end
        else
            try
                idnew = sub2ind(info.dim,id0(1),id0(2));
            catch
                ff=sprintf('please reset matrix index, maximum is [%d %d] \n',info.dim);
                error(ff)
            end
        end
        idnewf = find(idf==idnew, 1);
        fprintf('%d\n',idnewf)
        if isempty(idnewf)
            aa = sfprintf('matrix index error! data dimension is %dx%d, input ID is %d (=(%d-1)*%d+%d)\n',info.dim,idnew,id(2),info.dim(1),id(1));
            error(aa);
        end
    end
end
if flag_datatype0~=flag_datatype
    aa=sprintf('different data type! please check data information pannel\n;%d\n %d\n',flag_datatype0,flag_datatype);
    error(aa)
end
if ~strcmp(Datafold0,Datafold)
    bb=sprintf('different Data folder! please check data information pannel\n;%s\n %s\n',Datafold0,Datafold);
    error(bb)
end
if flag_datatype==1 && ~strcmp(ROI_mask0,ROI_mask)
    cc=sprintf('different Data folder! please check data information pannel\n;%s\n %s\n',ROI_mask0,ROI_mask);
    error(cc)
end

if ~isempty(setdiff(subid0,subid))
    dd = sprintf('Click Run then plot \n');
    error(dd)
end
wgr_plot_powerspectral(y,Fs,flag_datatype,matfull0,subid0,SubFold,idnewf)

function wgr_plot_powerspectral(y,Fs,flag_datatype,matfull,subid,SubFold,idnewf)
% t = 1:1000; % y = 3*sin(2*pi*0.05*t) + 2*randn(size(t));% TR = 1;% Fs = 1/TR; %sampling frequency
for isub=subid
    y0 = y{isub};
    if ~isempty(y0)
        y0 =y0(:,1);
        L = size(y0,1); %length of signal
        NFFTa(isub) = 2^nextpow2(L); % Next power of 2 from length of y
%         fprintf('NFFT=%d\n',NFFT)
    end
end
NFFT= max(NFFTa);
f = Fs/2*linspace(0,1,NFFT/2+1);

kk=1;
for isub=subid
    y0 = y{isub};
    if ~isempty(y0)
        y0 =y0(:,idnewf);       
        Y(:,kk) = fft(y0,NFFT)/L;
        lege{kk,1} = strrep(SubFold(isub).name,'_','-'); kk=kk+1;
    end
end
titname = sprintf('ROI-wise (#%d subjects)', nnz(subid));
wgr_plot_fft(titname,f,Y,NFFT,lege,[]);
        
if flag_datatype==2&&matfull==1
    kk=1;
    for isub=subid
        y0 = y{isub};
        if ~isempty(y0)
            tmp = fft(y0,NFFT)/L;
            Y(:,kk) = dynamicBC_nanmean(tmp,2); 
            err(:,kk)= std(2*abs(tmp(1:NFFT/2+1,:)),[],2);
            lege{kk,1} = strrep(SubFold(isub).name,'_','-'); kk=kk+1;
        end
    end
    titname = sprintf('Full matrix -- Mean spectrum (#%d subjects)', nnz(subid));
    wgr_plot_fft(titname,f,Y,NFFT,lege,err);
end

function []=wgr_select_dir_rewrite(varargin)
D = varargin{3};
str = varargin{4};
flag = get(D.tx_mask,'val');
if strcmp(str,'mask')
    if flag==1
        files = spm_select(1,'image','maskfile');
    end
else
    files = spm_select(1,'dir','file folder');
end

if strcmp(str ,'mask') && flag==1
    set(D.ed_mask,'string',files);
elseif strcmp(str ,'data')
    set(D.ed_data,'string',files);
elseif strcmp(str,'output')
    set(D.ed_outdir,'string',files);
end

function []=wgr_data_type(varargin)
D = varargin{3};  % Get structure.
flag = get(D.data_type,'val');
if flag==1
    set([D.tx_mask,D.ed_mask],'vis','on');   
    set([D.ed_varname,D.tx_varname],'vis','off');
    set([D.ed_ID,D.tx_ID],'vis','off');
    set(D.matrix_select,'vis','off');
    set(D.FCGC_select,'vis','on');
elseif flag==2
   set([D.tx_mask,D.ed_mask],'vis','off');
   set([D.ed_varname,D.tx_varname],'vis','on');
   set([D.ed_ID,D.tx_ID],'vis','on');
   set(D.matrix_select,'vis','on');
   set(D.FCGC_select,'vis','off');
else
    set([D.ed_varname,D.tx_varname],'vis','off');
    set([D.tx_mask,D.ed_mask],'vis','off');   
    set([D.ed_ID,D.tx_ID],'vis','off');
    set(D.matrix_select,'vis','off');
    set(D.FCGC_select,'vis','off');
end

function DynamicBC_Reslice(PI,PO,hld,TargetSpace)
%   PI - input filename
%   PO - output filename
%   hld - interpolation method. 0: Nearest Neighbour. 1: Trilinear.
%   TargetSpace - Define the target space. 'ImageItself': defined by the image itself (corresponds  to the new voxel size); 'XXX.img': defined by a target image 'XXX.img' (the NewVoxSize parameter will be discarded in such a case).

headIN = spm_vol(TargetSpace) ;
headIN = headIN(1);
mat=headIN.mat;
dim=headIN.dim;

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

function []=wgr_sub_select_call(varargin)
D = varargin{3};  % Get structure.
FP = get(D.sub_select,'val');
Datafold = get(D.ed_data,'string');
if ~exist(Datafold,'dir')
    return
end
SubFold = dir(Datafold); SubFold(1:2)=[];
for i=1:length(SubFold)
    str{i,1} = SubFold(i).name;
end
if FP==2
    sub.figselect = figure('Name','Select Subjects (Ctrl / Shift + Click)',...       
        'units','normalized',...      
        'menubar','none',...       
        'numbertitle','off',...      
        'unit','normalized',...
        'color',[0.95 0.95 0.95],...
        'position',[0.25 0.2 0.25 0.3]);
    sub.currentfolder = uicontrol('Parent',sub.figselect,'Style','pushbutton','units','normalized','fontunits', 'normalized', 'fontsize',0.5,...
       'String',Datafold, 'Position',[0 0.9 1 0.08],'enable','off');%,'CallBack','close(figselect)'
    sub.sublist = uicontrol('Parent',sub.figselect,'Style','Listbox','units','normalized','fontunits', 'normalized', 'fontsize',0.12,'BackgroundColor',[1 1 1],'Max',8, ...
       'String',str, 'Position',[0.05 0.05 0.65 0.8],'Tag','FileListLb','Value',[]);
    sub.okbutton = uicontrol('Parent',sub.figselect,'Style','pushbutton','units','normalized','fontunits', 'normalized', 'fontsize',0.5,...
       'String','Done', 'Position',[0.75 0.43 0.2 0.14]);%,'CallBack','close(figselect)'

    set(sub.okbutton,'callback',{@wgr_select_done_call,sub,D});
else
    subid = num2str([1:length(SubFold)]);
    set(D.sub_select,'TooltipString',subid) ;
end

function [subid]=wgr_select_done_call(varargin)
sub = varargin{3};  
D = varargin{4}; 
subid = num2str(get(sub.sublist,'val'));
if isempty(subid)
    Datafold = get(D.ed_data,'string');
    SubFold = dir(Datafold); SubFold(1:2)=[];
    subid = num2str([1:length(SubFold)]);
    set(D.sub_select,'TooltipString',subid) ;
end
set(D.sub_select,'TooltipString',subid) ;
close(sub.figselect)


function wgr_plot_fft(titname,f,Y,NFFT,lege,err)
fig= figure('Name',titname,...       
    'units','normalized',...      
    'numbertitle','off',...      
    'unit','normalized',...
    'color',[0.95 0.95 0.95],...
    'position',[0.25 0.2 0.5 0.3]);
% Plot single-sided amplitude spectrum.
axes1 = axes('Parent',fig,...
    'Position',[0.1 0.2 0.85 0.7],'fontunits','norm','fontsize',0.07);
box(axes1,'on');
hold(axes1,'all');
if isempty(err)
    plot(f,2*abs(Y(1:NFFT/2+1,:))) 
else
    colo = jet(size(Y,2));
    hold all
    for i=1:size(Y,2)
        errorbar(f,2*abs(Y(1:NFFT/2+1,i)),err(:,i),'-',...
            'Color',colo(i,:));
    end
    xlim(f([1 end]))
    hold off
end
legend(lege)
title('Single-Sided Amplitude Spectrum of y(t)','fontunits','norm','fontsize',0.07)
xlabel('Frequency (Hz)','fontunits','norm','fontsize',0.07)
ylabel('|Y(f)|','fontunits','norm','fontsize',0.07);