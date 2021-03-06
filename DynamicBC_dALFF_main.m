function [] = DynamicBC_dALFF_main
D.fig = figure('Name','dynamic ALFF analysis',...       
    'units','normalized',...      
    'menubar','none',...       
    'numbertitle','off',...      
    'unit','normalized',...
    'color',[0.95 0.95 0.95],...
    'position',[0.25 0.2 0.5 0.4]);
movegui(D.fig,'center');      

ImgInfopos = [0.1 0.5 0.38 0.4];
D.ImgInfo = uibuttongroup('units','norm',...
    'Title','Information for alff Process','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
    'pos',ImgInfopos);

Dyinfopos = [0.52 0.5 0.38 0.4];
D.DyInfo = uibuttongroup('units','norm',...
    'Title','Dynamic Parameter','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
    'pos',Dyinfopos);
              
IOpos = [0.1 0.1 0.6 0.35];
D.IO = uibuttongroup('units','norm',...
    'Title','In/Out Information','fontunits', 'normalized', 'fontsize',0.1,'foregroundcolor',[.1 .1 .1],...
    'pos',IOpos);

D.Runbut = uicontrol('Parent',D.fig,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',[0.75 0.3 0.15 0.1],...
                    'string','Run',...       
                    'fontunits', 'normalized', 'fontsize',0.6,...
                    'fontweight','bold',...
                    'horizontalalign','center');

D.Returnbut = uicontrol('Parent',D.fig,...
                    'style','pushbutton',...
                    'units','norm',...
                    'position',[0.75 0.15 0.15 0.1],...
                    'string','Return',...       
                    'fontunits', 'normalized', 'fontsize',0.6,...
                    'fontweight','bold',...
                    'horizontalalign','center');
%% MASK
D.tx_mask = uicontrol('Parent',D.IO,...
                'style','text',...
                'units','norm',...
                'position',[0.02 0.66 0.22 0.3],...
                'string','Mask files',...
                'backgroundc',get(D.IO,'backgroundc'),...
                'fontunits', 'normalized', 'fontsize',0.5,...
                'fontweight','bold',...
                'horizontalalign','left',...
                'foregroundcolor',[0 0 0.8]);
D.ed_mask = uicontrol('Parent',D.IO,...
                        'style','edit',...
                        'unit','norm',...
                        'position',[0.25 0.66 0.6 0.3],...
                        'fontunits', 'normalized', 'fontsize',0.6,...
                        'string','not selected');  
D.ed_maskselect = uicontrol('Parent',D.IO,...
                            'style','pushbutton',...
                            'units','norm',...
                            'position',[0.88 0.66 0.09 0.3],...
                            'string','...',...
                            'TooltipString','Click me!',...
                            'fontunits', 'normalized', 'fontsize',0.8);  

%% IN
D.tx_IN = uicontrol('Parent',D.IO,...
                'style','text',...
                'units','norm',...
                'position',[0.02 0.34 0.22 0.3],...
                'string','Input',...
                'backgroundc',get(D.IO,'backgroundc'),...
                'fontunits', 'normalized', 'fontsize',0.5,...
                'fontweight','bold',...
                'horizontalalign','left',...
                'foregroundcolor',[0 0 0.8]);
D.ed_IN = uicontrol('Parent',D.IO,...
                        'style','edit',...
                        'unit','norm',...
                        'position',[0.25 0.34 0.6 0.3],...
                        'fontunits', 'normalized', 'fontsize',0.6,...
                        'string','not selected');  
D.ed_INselect = uicontrol('Parent',D.IO,...
                            'style','pushbutton',...
                            'units','norm',...
                            'position',[0.88 0.34 0.09 0.3],...
                            'string','...',...
                            'TooltipString','Click me!',...
                            'fontunits', 'normalized', 'fontsize',0.8);  
%% OUT
D.tx_OUT = uicontrol('Parent',D.IO,...
                'style','text',...
                'units','norm',...
                'position',[0.02 0.02 0.22 0.3],...
                'string','Output',...
                'backgroundc',get(D.IO,'backgroundc'),...
                'fontunits', 'normalized', 'fontsize',0.5,...
                'fontweight','bold',...
                'horizontalalign','left',...
                'foregroundcolor',[0 0 0.8]);
D.ed_OUT = uicontrol('Parent',D.IO,...
                        'style','edit',...
                        'unit','norm',...
                        'position',[0.25 0.02 0.6 0.3],...
                        'fontunits', 'normalized', 'fontsize',0.6,...
                        'string','not selected');  
D.ed_OUTselect = uicontrol('Parent',D.IO,...
                            'style','pushbutton',...
                            'units','norm',...
                            'position',[0.88 0.02 0.09 0.3],...
                            'string','...',...
                            'TooltipString','Click me!',...
                            'fontunits', 'normalized', 'fontsize',0.8);  
 %% TR for filter                 
D.TRtxt = uicontrol('Parent',D.ImgInfo,...
    'style','text',...
    'units','norm',...
    'pos',[0.1,0.66,0.55,0.3],...
    'string','TR(s)',...
    'fontunits', 'normalized', 'fontsize',0.45);
D.TRedit = uicontrol('Parent',D.ImgInfo,...
    'style','edit',...
    'units','norm',...
    'pos',[0.65,0.66,0.25,0.3],...
    'fontunits', 'normalized', 'fontsize',0.45,...
    'string','2');  
%% Lowcut
D.Lowtxt = uicontrol('Parent',D.ImgInfo,...
    'style','text',...
    'units','norm',...
    'pos',[0.1,0.34,0.55,0.3],...
    'string','LowCut(Hz)',...
    'fontunits', 'normalized', 'fontsize',0.45);
D.Lowedit = uicontrol('Parent',D.ImgInfo,...
    'style','edit',...
    'units','norm',...
    'pos',[0.65,0.34,0.25,0.3],...
    'fontunits', 'normalized', 'fontsize',0.45,...
    'string','0.01');   
%% Highcut
D.Hightxt = uicontrol('Parent',D.ImgInfo,...
    'style','text',...
    'units','norm',...
    'pos',[0.1,0.02,0.55,0.3],...
    'string','HighCut(Hz)',...
    'fontunits', 'normalized', 'fontsize',0.45);
D.Highedit = uicontrol('Parent',D.ImgInfo,...
    'style','edit',...
    'units','norm',...
    'pos',[0.65,0.02,0.25,0.3],...
    'fontunits', 'normalized', 'fontsize',0.45,...
    'string','0.08');
%% window length
D.Winlen = uicontrol('Parent',D.DyInfo,...
    'style','text',...
    'units','norm',...
    'pos',[0.1 0.55 0.50 0.4],...
    'string','WinLen',...
    'fontunits','normalized', 'fontsize',0.45);
D.Winlened = uicontrol('Parent',D.DyInfo,...
    'style','edit',...
    'units','norm',...
    'pos',[0.65 0.55 0.25 0.4],...
    'string','50',...
    'fontunits','normalized', 'fontsize',0.45);
%% ovelap
D.overlap = uicontrol('Parent',D.DyInfo,...
    'style','text',...
    'units','norm',...
    'pos',[0.1 0.05 0.50 0.4],...
    'string','Overlap',...
    'fontunits','normalized', 'fontsize',0.45);
D.overlaped = uicontrol('Parent',D.DyInfo,...
    'style','edit',...
    'units','norm',...
    'pos',[0.65 0.05 0.25 0.4],...
    'string','0.6',...
    'fontunits','normalized', 'fontsize',0.45);

set(D.ed_maskselect,'callback',{@dbc_selmask,D});
set(D.ed_INselect,'callback',{@dbc_selInput,D});
set(D.ed_OUTselect,'callback',{@dbc_selOutput,D});
set(D.Runbut,'callback',{@dbc_dalffrun,D});
set(D.Returnbut,'callback',{@dbc_return,D});

end

function dbc_selmask(varargin)
D = varargin{3};
files = spm_select(1,'image','maskfile');
set(D.ed_mask,'string',files(1:end-2))
end

function dbc_selInput(varargin)
D = varargin{3};
files = spm_select(1,'dir','Input with subfolds');
set(D.ed_IN,'string',files)
end

function dbc_selOutput(varargin)
D = varargin{3};
files = spm_select(1,'dir','Output dirs');
set(D.ed_OUT,'string',files)
end

function dbc_dalffrun(varargin)
D = varargin{3};
TRs = str2num(get(D.TRedit,'string'));
Highcut = str2num(get(D.Highedit,'string'));
Lowcut = str2num(get(D.Lowedit,'string'));
Winlen = str2num(get(D.Winlened,'string'));
Overlap = str2num(get(D.overlaped,'string'));

save_info.TR = TRs;
save_info.highcut = Highcut;
save_info.lowcut = Lowcut;
save_info.slw_alignment = 1;

maskfiles = get(D.ed_mask,'string');
INdirs = get(D.ed_IN,'string');
Outdiro = get(D.ed_OUT,'string');
Outdir1d = fullfile(Outdiro,'dALFFmap');
Outdir1m = fullfile(Outdiro,'zdALFFmap');
Outdir1z = fullfile(Outdiro,'mdALFFmap');
Outdir2 = fullfile(Outdiro,'variance_dALFFmap');
subnames = dir(INdirs);
[vm datm] = Dynamic_read_dir_NIFTI(maskfiles);
indexs = find(datm);
% save test
for isub = 1:length(subnames)-2
    mkdir(fullfile(Outdir1d,subnames(isub+2).name));
    mkdir(fullfile(Outdir1m,subnames(isub+2).name));
    mkdir(fullfile(Outdir1z,subnames(isub+2).name));
    mkdir(fullfile(Outdir2,subnames(isub+2).name));
    
    [vo dato] = Dynamic_read_dir_NIFTI(fullfile(INdirs,subnames(isub+2).name));
    data = dato(indexs,:)';
    [dALFF] = DynamicBC_dALFF(data,Winlen,Overlap,save_info);
%     save(subnames(isub+2).name,'dALFF')
    for i = 1:length(dALFF)
        DAT0 = zeros(vo(1).dim);
        zDAT0 = DAT0;
        mDAT0 = DAT0;
        DAT0(indexs) = dALFF{i};
        zdALFF(i,:) = (dALFF{i}-mean(dALFF{i}))/std(dALFF{i});
        mdALFF(i,:) = dALFF{i}/mean(dALFF{i});
        zDAT0(indexs) = (dALFF{i}-mean(dALFF{i}))/std(dALFF{i});
        mDAT0(indexs) = dALFF{i}/mean(dALFF{i});
        
        if i<10
            fname = fullfile(fullfile(Outdir1d,subnames(isub+2).name),['dALFF_000',num2str(i),'.nii']);
            zfname = fullfile(fullfile(Outdir1z,subnames(isub+2).name),['zdALFF_000',num2str(i),'.nii']);
            mfname = fullfile(fullfile(Outdir1m,subnames(isub+2).name),['mdALFF_000',num2str(i),'.nii']);
        elseif i<100
            fname = fullfile(fullfile(Outdir1d,subnames(isub+2).name),['dALFF_00',num2str(i),'.nii']);
            zfname = fullfile(fullfile(Outdir1z,subnames(isub+2).name),['zdALFF_00',num2str(i),'.nii']);
            mfname = fullfile(fullfile(Outdir1m,subnames(isub+2).name),['mdALFF_00',num2str(i),'.nii']);
        elseif i<1000
            fname = fullfile(fullfile(Outdir1d,subnames(isub+2).name),['dALFF_0',num2str(i),'.nii']);
            zfname = fullfile(fullfile(Outdir1z,subnames(isub+2).name),['zdALFF_0',num2str(i),'.nii']);
            mfname = fullfile(fullfile(Outdir1m,subnames(isub+2).name),['mdALFF_0',num2str(i),'.nii']);
        elseif i<10000
            fname = fullfile(fullfile(Outdir1d,subnames(isub+2).name),['dALFF_',num2str(i),'.nii']);
            zfname = fullfile(fullfile(Outdir1z,subnames(isub+2).name),['zdALFF_',num2str(i),'.nii']);
            mfname = fullfile(fullfile(Outdir1m,subnames(isub+2).name),['mdALFF_',num2str(i),'.nii']);
        else
            error('too many output');
        end
        DynamicBC_write_NIFTI(DAT0,vo,fname);
        DynamicBC_write_NIFTI(zDAT0,vo,zfname);
        DynamicBC_write_NIFTI(mDAT0,vo,mfname);
    end
    for i = 1:length(dALFF)
        DALFF(i,:) = dALFF{i};
    end
    for i = 1:size(DALFF,2)
        stddalff(1,i) = std(DALFF(:,i));
        stdzdalff(1,i) = std(zdALFF(:,i));
        stdmdalff(1,i) = std(mdALFF(:,i));
    end
    fname = fullfile(fullfile(Outdir2,subnames(isub+2).name),['var_dALFF.nii']);
    zfname = fullfile(fullfile(Outdir2,subnames(isub+2).name),['var_zdALFF.nii']);
    mfname = fullfile(fullfile(Outdir2,subnames(isub+2).name),['var_mdALFF.nii']);
    
    DAT0 = zeros(vo(1).dim);
    DAT0(indexs) = stddalff;
    DynamicBC_write_NIFTI(DAT0,vo,fname);
    
    DAT0 = zeros(vo(1).dim);
    DAT0(indexs) = stdzdalff;
    DynamicBC_write_NIFTI(DAT0,vo,zfname);
    
    DAT0 = zeros(vo(1).dim);
    DAT0(indexs) = stdmdalff;
    DynamicBC_write_NIFTI(DAT0,vo,mfname);
    disp([subnames(isub+2).name,' finished'])
end

disp('dALFF all finished')
end

function dbc_return(varargin)
D = varargin{3};
close(D.fig);
end