function DynamicBC_Cluster()
%version 3.0 2018.08.02
D.fig = figure('Name','K-means Clustering (v3.0)',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'unit','normalized',...
    'color',[0.95 0.95 0.95],...
    'position',[0.25 0.2 0.5 0.2]);
movegui(D.fig,'center');

pos_Data = [0.01 0.05 0.65 0.9];
pos_Clustering = [0.68 0.05 0.31 0.9];
             
D.data_info = uibuttongroup('units','norm',...
    'Title','Data Information','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
                     'pos',pos_Data);   

pos_data_type = [0.05 0.71 0.49 0.2];
pos_conn_type = [0.56 0.71 0.39 0.2];
pos_data_mask_str  =  [0.05 0.49 0.35 0.2];
pos_data_mask_str_input  =  [0.4 0.49 0.55 0.2];
pos_data_select_str = [0.05 0.27 0.35 0.2];
pos_data_select = [0.4 0.27 0.55 0.2];
pos_out_str = [0.05 0.05 0.35 0.2];
pos_out_select = [0.4 0.05 0.55 0.2];

%% FC/GC selection
D.FG = uicontrol('Parent',D.data_info,...
                'style','popup',...
                'units','norm',...
                'position',pos_conn_type,...
                'string',{'FC/dALFF','GC'},...               
                'fontunits', 'normalized', 'fontsize',0.5,...
                'foregroundcolor',[.1 .1 .1],...
                'value',1); 
            
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
% NIFTI image
D.tx_mask = uicontrol('Parent',D.data_info,...
                'style','pushbutton',...
                'units','norm',...
                'position',pos_data_mask_str,...
                'string','Mask(NIFTI):',...
                'backgroundc',get(D.data_info,'backgroundc'),...                
                'visible','off',...
                'TooltipString','Click me to select NIFTI image',...
                'fontunits', 'normalized', 'fontsize',0.5,...              
                'horizontalalign','left');    
                
D.ed_mask = uicontrol('Parent',D.data_info,...
                        'style','edit',...
                        'unit','norm',...
                        'position',pos_data_mask_str_input,...
                        'fontunits', 'normalized', 'fontsize',0.5,...                        
                        'visible','off',...
                        'string','',...
                        'horizontalalign','right');  
              
%% Matrix                        
D.tx_varname = uicontrol('Parent',D.data_info,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_data_mask_str,...
                    'string','Variable:',...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...                   
                    'value',0,...                    
                    'visible','off',...                   
                    'TooltipString','e.g. FCM.Matrix');                 
D.ed_varname = uicontrol('Parent',D.data_info,...
                    'style','edit',...
                    'units','norm',...
                    'position',pos_data_mask_str_input,...
                    'string','FCM.Matrix',...                    
                    'backgroundc',get(D.data_info,'backgroundc'),...
                    'fontunits', 'normalized', 'fontsize',0.5,...               
                    'horizontalalign','center',...                   
                    'visible','off',...
                    'TooltipString','FCM.Matrix or GCM.Matrix');        
                        
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
              
%% Cluster methods selection
pos_CMe = [0.05 0.75 0.9 0.2] ;
pos_clusize_str = [0.05 0.5 0.6 0.2] ;
pos_clusize_input = [0.65 0.5 0.3 0.2] ;
pos_Run = [0.05 0.27 0.9 0.2] ;
pos_Plot = [0.05 0.05 0.9 0.2] ;

D.Clustering = uibuttongroup('Parent',D.fig,...
                     'units','norm',...
                     'Title','Clustering','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
                     'pos',pos_Clustering);   
D.CMe = uicontrol('Parent',D.Clustering,...
                'style','popup',...
                'units','norm',...
                'position',pos_CMe,...
                'string',{'sqeuclidean','correlation','cityblock','Distance Measures'},...               
                'fontunits', 'normalized', 'fontsize',0.5,...                
                'foregroundcolor',[.1 .1 .1],...
                'TooltipString','Choose one of them (1~3)',...
                'value',4); 

D.tx_fs = uicontrol('Parent',D.Clustering,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',pos_clusize_str,...
                    'string','Cluster #:',...                                      
                    'TooltipString','Click me to estimate number of clusters, all results will be saved (fixed/estimated)',...
                    'fontunits', 'normalized', 'fontsize',0.5,...                   
                    'horizontalalign','center');  
D.ed_fs = uicontrol('Parent',D.Clustering,...
                    'style','edit',...
                    'units','norm',...
                    'position',pos_clusize_input,...
                    'string','6',...                   
                    'fontunits', 'normalized', 'fontsize',0.5,...              
                    'horizontalalign','center',...                   
                    'visible','on',...
                    'TooltipString','Fixed number of clusters');   

D.pb_run = uicontrol(D.Clustering,...
                  'style','pushbutton',...
                  'units','norm',...
                  'position',pos_Run,...
                  'string','RUN',...                
                  'TooltipString','Click me!',...
                  'fontunits', 'normalized', 'fontsize',0.6,'val',1); 
D.pb_plot = uicontrol(D.Clustering,...
                  'style','pushbutton',...
                  'units','norm',...
                  'position',pos_Plot,...
                  'string','Results plot',...
                  'TooltipString','Click me After RUN!',...
                  'fontunits', 'normalized', 'fontsize',0.6);     
              
set(D.FG,'callback',{@wgr_FCGC_select,D});
set(D.tx_mask,'callback',{@wgr_select_dir_rewrite,D,'mask',''});
set(D.data_tx,'callback',{@wgr_select_dir_rewrite,D,'data',''});
set(D.out_str,'callback',{@wgr_select_dir_rewrite,D,'output',''});
set(D.data_type,'callback',{@wgr_data_type,D});
set(D.tx_fs,'callback',{@wgr_estimatClu_call,D});
set(D.pb_run,'callback',{@wgr_run_call,D});
set(D.pb_plot,'callback',{@wgr_plot_call,D});

end

function []=wgr_FCGC_select(varargin)
D = varargin{3};  % Get structure.
FG = get(D.FG,'val');
if FG==1
    set(D.ed_varname,'string','FCM.Matrix');
else
    set(D.ed_varname,'string','GCM.Matrix');
end
end

function wgr_estimatClu_call(varargin)
D = varargin{3};  % Get structure.
sty = get(D.tx_fs,'style');
switch sty;
        case 'pushbutton'
            set(D.tx_fs,'style','popup');
            set(D.tx_fs,'string',{'Fixed','Estimate'});
end

end

function []=wgr_run_call(varargin)
D = varargin{3};  % Get structure.
flag_datatype = get(D.data_type,'val');
Datafold = get(D.ed_data,'string');
SubFold = dir(Datafold);
outputd = get(D.ed_outdir,'string');
k = str2num(get(D.ed_fs,'string'));  % cluster number
k = unique(k);
Distid = get(D.CMe,'val');
Distmethod0 = get(D.CMe,'string');
if Distid==length(Distmethod0)
    Distmethod = Distmethod0{1};
else
    Distmethod = Distmethod0{Distid};
end

FG = get(D.FG,'val');
flag_estimate = get(D.tx_fs,'val')-1;
if FG==1
    disp('FC/dALFF Clustering')
elseif FG==2&&flag_datatype==1
    disp('GC Maps Clustering (IN/OUT)')
else
    disp('GC Matrices Clustering')
end

if ~exist(outputd,'dir')
    mkdir(outputd)
end
if flag_datatype==1 %NIFTI
    maskfile = get(D.ed_mask,'string');
    if FG==1
        DynamicBC_clustermaps(k,outputd,maskfile,Datafold,Distmethod,'',flag_estimate);
    else
        DynamicBC_clustermaps(k,outputd,maskfile,Datafold,Distmethod,'*_IN_',flag_estimate);
        DynamicBC_clustermaps(k,outputd,maskfile,Datafold,Distmethod,'*_OUT_',flag_estimate);
    end
else
    matname = get(D.ed_varname,'string');
    DynamicBC_clustermatrix(k,outputd,matname,Datafold,Distmethod,flag_estimate);
end

fprintf('Done!  =_= \n')
end

function []=wgr_plot_call(varargin)
D = varargin{3};  % Get structure.
fpath = get(D.ed_outdir,'string');
matfile0 = spm_select('FPListRec',fpath,'^cluster_index.*\.mat$')
ll = length('cluster_index');
for imat=1:size(matfile0,1)
    matfile = strcat(matfile0(imat,:));
    [fpath0,name,ext] = fileparts(matfile);
    if length(name)>ll
        strio = [name(ll+2:end)];
    else
        strio = '';
    end
    fprintf('load %s\n',matfile)
    load(matfile);
    fig = figure('Name','State Transition',...       
        'units','normalized',...      
        'unit','normalized',...
        'color',[0.95 0.95 0.95],...
        'position',[0.1 0.2 0.5 0.3]);

        axes('parent',fig,'pos',[0.1,0.15,0.8,0.75]);
        IDX_subj0 = [];
        nsub = length(IDX_subj);
        for ii=1:nsub
            IDX_subj0(ii,1:length(IDX_subj{ii})) = IDX_subj{ii};
        end
        imagesc(IDX_subj0)
        if ~isempty(strio)
            title(['State Transition (K=',num2str(K),', ',strio,')'])
        else
            title(['State Transition (K=',num2str(K),')'])
        end
        xlabel('Time')
        ylabel('Subjects')
        colorbar('eastoutside');
        TM0 = 0;
        for ii=1:nsub
            [F{ii}, TM{ii}, MDT{ii}, NT{ii}] = icatb_dfnc_statevector_stats(IDX_subj{ii}, K);
            TM0 = TM{ii}+TM0;
        end
        TM0 = TM0./nsub;
        save(fullfile(fpath0,['cluster_stats',strio,'.mat']),'F','TM','MDT','NT')
        F0 = cell2mat(F');
        MDT0 = cell2mat(MDT');
        NT0 = cell2mat(NT);
        
    fig2 = figure('Name',['State Transition (K=',num2str(K),',',strio,')'],...       
        'units','normalized',...      
        'unit','normalized',...
        'color',[0.95 0.95 0.95],...
        'position',[0.1 0.2 0.6 0.3]);
        h1= subplot(1,19,1:5); 
        m = mean(F0);
        s = std(F0)/sqrt(nsub);
        errorbar(1:K,m,s,'-s','MarkerSize',10,...
            'MarkerEdgeColor','red','MarkerFaceColor','red')
        axis square
        xlim([0.8 K+0.2])
        ms = [m-s m+s]; msd = range(ms)/20;
        ylim([min(ms)-msd max(ms)+msd])
        set(h1,'xtick',1:K)
        xlabel('State (Cluster index)'); 
        ylabel('Frequency (mean/SE)')
        
        h2 = subplot(1,19,7:11);
        m = mean(MDT0);
        s = std(MDT0)/sqrt(nsub);
        errorbar(1:K,m,s,'-s','MarkerSize',10,...
            'MarkerEdgeColor','red','MarkerFaceColor','red')
        axis square
        set(h2,'xtick',1:K)
        xlim([0.8 K+0.2])
        ms = [m-s m+s]; msd = range(ms)/20;
        ylim([min(ms)-msd max(ms)+msd])
        xlabel('State (Cluster index)'); 
        ylabel('Mean dwell time (mean/SE)')
        
        h3 = subplot(1,19,13:19); 
        imagesc(TM0);  axis square
        set(h3,'xtick',1:K,'ytick',1:K)
        xlabel('State at T'); 
        ylabel('State at T-1')
        hC = colorbar('westoutside');
        ylabel(hC,'Probability','FontSize',16)
        caxis([0 1])
        NIstr = sprintf('Transitions #: %3.1f/%3.1f(M/SD)',mean(NT0),std(NT0));%         NIstr = sprintf('Number of transitions: %3.1f/%3.1f(M/SD)',mean(NT0),std(NT0));
        title(NIstr)
        
    a=spm_select('FPListrec',fpath0,'^Centroid_.*\.nii$'); % do not show for GC
    if ~isempty(a)    
        spm_check_registration(a);
        colormap(jet);
    end
%% plot matrix cluster...
    if exist(fullfile(fpath0,'Centroid_1.mat'),'file')
        matname={};
        for i=1:K
            matname{i,1} = fullfile(fpath0,['Centroid_',num2str(i),'.mat']);
        end
        for i=1:K
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
        nobs= length(IDXall);
        figure('Name',sprintf('Centroid (total #%d maps)',nobs),'Position',[100 scrsz(4)/3 scrsz(3)*4/5 scrsz(4)/4])

        for i=1:K
            clear DAT
            load(matname{i})
            DAT(1:size(DAT,1)+1:end)=0;
            subplot(1,K,i);imagesc(DAT);
            axis square;axis off
            set(gca, 'CLim', [mn mx]);
            freqall = nnz(IDXall==i);
            strfreq = sprintf('%d (#%d maps,%2.2f%%)',i, freqall,freqall/nobs*100);
            title(strfreq)
        end
        colormap(jet);
    end
end
end

function []=wgr_select_dir_rewrite(varargin)
D = varargin{3};
str = varargin{4};
if strcmp(str,'mask')
    files = spm_select(1,'image','maskfile');
else
    files = spm_select(1,'dir','file folder');
end

if strcmp(str ,'mask')
    set(D.ed_mask,'string',files);
elseif strcmp(str ,'data')
    set(D.ed_data,'string',files);
elseif strcmp(str,'output')
    set(D.ed_outdir,'string',files);
end
end

function []=wgr_data_type(varargin)
D = varargin{3};  % Get structure.
flag = get(D.data_type,'val');
if flag==1
    set([D.tx_mask,D.ed_mask],'vis','on');   
    set([D.ed_varname,D.tx_varname],'vis','off');
elseif flag==2
   set([D.tx_mask,D.ed_mask],'vis','off');
   set([D.ed_varname,D.tx_varname],'vis','on');
else
    set([D.ed_varname,D.tx_varname],'vis','off');
    set([D.tx_mask,D.ed_mask],'vis','off');   
end
end

function [F, TM, MDT, NT] = icatb_dfnc_statevector_stats(a, k)

Nwin = length(a);

%% Fraction of time spent in each state
F = zeros(1,k);
for jj = 1:k
    F(jj) = (sum(a == jj))/Nwin;
end

%% Number of Transitions
NT = sum(abs(diff(a)) > 0);

%% Mean dwell time in each state
MDT = zeros(1,k);
for jj = 1:k
    start_t = find(diff(a==jj) == 1);
    end_t = find(diff(a==jj) == -1);
    if a(1)==jj
        start_t = [0; start_t];
    end
    if a(end) == jj
        end_t = [end_t; Nwin];
    end
    MDT(jj) = mean(end_t-start_t);
    if isempty(end_t) & isempty(start_t)
        MDT(jj) = 0;
    end
end

%% Full Transition Matrix
TM = zeros(k,k);
for t = 2:Nwin
    TM(a(t-1),a(t)) =  TM(a(t-1),a(t)) + 1;
end

for jj = 1:k
    if sum(TM(jj,:)>0)
        TM(jj,:) = TM(jj,:)/sum(a(1:Nwin-1) == jj);
    else
        TM(jj,jj) = 1;
    end
end

end