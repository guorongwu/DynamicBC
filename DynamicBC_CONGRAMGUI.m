function DynamicBC_CONGRAMGUI
% Hsize = get(0,'screensize');
[pat,nam,ext] = fileparts(which('DynamicBC_CONGRAMGUI.m'));
load(fullfile(pat,'COLORMODULE.mat'));
if ~isempty(dir(fullfile(pat,'COLTEMP.mat')))
    delete(fullfile(pat,'COLTEMP.mat'));
end
DBCCGG.fig = figure('units','norm','pos',[0.1,0.1,0.6,0.6],'name','Connectogram GUI v1.0','menu','none');
DBCCGG.IO = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.1,0.75,0.8,0.2]);
DBCCGG.Mod = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.1,0.65,0.8,0.1]);
DBCCGG.modax = axes('parent',DBCCGG.fig,'units','norm','pos',[0.1,0.55,0.8,0.1]);
Modc = 1:18;
image(Modc,'parent',DBCCGG.modax,'CDataMapping','scaled');colormap(DBCCGG.modax,COLORMODULE);axis(DBCCGG.modax,'off');
DBCCGG.nodeproperty = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.1,0.15,0.5,0.4],'title','Node Property Sel');
DBCCGG.Edgethr = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.6,0.35,0.3,0.2],'title','edge threshold');
DBCCGG.Degthr = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.6,0.25,0.3,0.1],'title','node threshold');
DBCCGG.Edgewidth = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.6,0.15,0.3,0.1],'title','edge width');
DBCCGG.Nodsize = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.1,0.05,0.2,0.1],'title','node size');
DBCCGG.Overlap = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.3,0.05,0.2,0.1],'title','Edge Above');
DBCCGG.Edgeshowtype = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.5,0.05,0.2,0.1],'title','Between Edge color');
DBCCGG.Opt = uibuttongroup('parent',DBCCGG.fig,'units','norm','pos',[0.7,0.05,0.2,0.1]);
DBCCGG.IO_inputtxt = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.1,0.55,0.1,0.4],'style','text','string','InputMatrix');
DBCCGG.IO_inputedi = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.2,0.55,0.6,0.4],'style','edit','string','NULL');
DBCCGG.IO_inputsel = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.8,0.55,0.1,0.4],'style','pushbutton','string','...');
DBCCGG.Mod_inputtxt = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.1,0.05,0.1,0.4],'style','text','string','InputModule');
DBCCGG.Mod_inputedi = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.2,0.05,0.6,0.4],'style','edit','string','NULL');
DBCCGG.Mod_inputsel = uicontrol('parent',DBCCGG.IO,'units','norm','pos',[0.8,0.05,0.1,0.4],'style','pushbutton','string','...');
DBCCGG.modc_txt = uicontrol('parent',DBCCGG.Mod,'units','norm','pos',[0.1,0.05,0.4,0.9],'style','text','string','Module number(max 18):');
DBCCGG.modc_txtnum = uicontrol('parent',DBCCGG.Mod,'units','norm','pos',[0.6,0.05,0.1,0.9],'style','text','string','18');
DBCCGG.modc_sel = uicontrol('parent',DBCCGG.Mod,'units','norm','pos',[0.7,0.05,0.2,0.9],'style','pushbutton','string','module color');
DBCCGG.nodeproperty_sel = uicontrol('parent',DBCCGG.nodeproperty,'units','norm','pos',[0.05,0.8,0.9,0.15],'style','pushbutton','string','select node property');
DBCCGG.nodeproperty_listbox = uicontrol('parent',DBCCGG.nodeproperty,'units','norm','pos',[0.05,0.05,0.45,0.7],'style','listbox');
DBCCGG.nodeproperty_listbox2 = uicontrol('parent',DBCCGG.nodeproperty,'units','norm','pos',[0.5,0.05,0.45,0.7],'style','listbox');
set(DBCCGG.nodeproperty_listbox2,'enable','off')
DBCCGG.Edgethr_txt = uicontrol('parent',DBCCGG.Edgethr,'units','norm','pos',[0.05,0.5,0.25,0.45],'style','text','string','threshold');
DBCCGG.Edgethr_edi = uicontrol('parent',DBCCGG.Edgethr,'units','norm','pos',[0.05,0.05,0.25,0.45],'style','edit','string','0.3');
DBCCGG.Edgethr_TYPE = uibuttongroup('parent',DBCCGG.Edgethr,'units','norm','pos',[1/3,0,1/3,1]);
DBCCGG.Edgethr_cost= uicontrol('parent',DBCCGG.Edgethr_TYPE,'units','norm','pos',[0.05,0.05,0.9,0.45],'style','rad','string','Cost','val',0);
DBCCGG.Edgethr_T = uicontrol('parent',DBCCGG.Edgethr_TYPE,'units','norm','pos',[0.05,0.5,0.9,0.45],'style','rad','string','T','val',1);
DBCCGG.Edgethr_AP = uibuttongroup('parent',DBCCGG.Edgethr,'units','norm','pos',[2/3,0,1/3,1]);
DBCCGG.Edgethr_abs = uicontrol('parent',DBCCGG.Edgethr_AP,'units','norm','pos',[0.05,0.05,0.9,0.45],'style','rad','string','ABS','val',0);
DBCCGG.Edgethr_pos = uicontrol('parent',DBCCGG.Edgethr_AP,'units','norm','pos',[0.05,0.5,0.9,0.45],'style','rad','string','POS','val',1);
DBCCGG.Degthr_hub = uicontrol('parent',DBCCGG.Degthr,'units','norm','pos',[0.05,0.05,0.3,0.9],'style','rad','string','>mean+std','val',1);
DBCCGG.Degthr_othertxt = uicontrol('parent',DBCCGG.Degthr,'units','norm','pos',[0.35,0.05,0.3,0.9],'style','rad','string','Threshold','val',0);
DBCCGG.Degthr_otheredi = uicontrol('parent',DBCCGG.Degthr,'units','norm','pos',[0.65,0.05,0.3,0.9],'style','edit','string','NULL');
set(DBCCGG.Degthr_otheredi,'enable','off');
DBCCGG.Edgewidth_equ = uicontrol('parent',DBCCGG.Edgewidth,'units','norm','pos',[0.05,0.05,0.3,0.9],'style','rad','string','equal','val',1);
DBCCGG.Edgewidth_cha = uicontrol('parent',DBCCGG.Edgewidth,'units','norm','pos',[0.35,0.05,0.3,0.9],'style','rad','string','raw','val',0);
DBCCGG.Edgewidth_chatxt = uicontrol('parent',DBCCGG.Edgewidth,'units','norm','pos',[0.65,0.05,0.2,0.9],'style','text','string','factor');
DBCCGG.Edgewidth_chaedi = uicontrol('parent',DBCCGG.Edgewidth,'units','norm','pos',[0.85,0.05,0.1,0.9],'style','edit','string','2');
DBCCGG.Nodsize_txt = uicontrol('parent',DBCCGG.Nodsize,'units','norm','pos',[0.05,0.05,0.45,0.9],'style','text','string','Node size');
DBCCGG.Nodsize_edi = uicontrol('parent',DBCCGG.Nodsize,'units','norm','pos',[0.5,0.05,0.45,0.9],'style','edit','string','8');
DBCCGG.Overlap_with = uicontrol('parent',DBCCGG.Overlap,'units','norm','pos',[0.05,0.05,0.45,0.9],'style','rad','string','Within','val',1);
DBCCGG.Overlap_bet = uicontrol('parent',DBCCGG.Overlap,'units','norm','pos',[0.5,0.05,0.45,0.9],'style','rad','string','Between','val',0);
DBCCGG.Edgeshowtype_gray = uicontrol('parent',DBCCGG.Edgeshowtype,'units','norm','pos',[0.05,0.05,0.45,0.9],'style','rad','string','gray','val',1);
DBCCGG.Edgeshowtype_comb = uicontrol('parent',DBCCGG.Edgeshowtype,'units','norm','pos',[0.5,0.05,0.45,0.9],'style','rad','string','combined','val',0);
DBCCGG.Opt_show = uicontrol('parent',DBCCGG.Opt,'units','norm','pos',[0.025,0.05,0.45,0.9],'style','pushbutton','string','Show');
DBCCGG.Opt_exit = uicontrol('parent',DBCCGG.Opt,'units','norm','pos',[0.525,0.05,0.45,0.9],'style','pushbutton','string','Exit');
set(DBCCGG.IO_inputsel,'callback',{@DBCCGGinputmat,DBCCGG});
set(DBCCGG.Mod_inputsel,'callback',{@DBCCGGinputmod,DBCCGG});
set(DBCCGG.modc_sel,'callback',{@DBCCGGmodcol,DBCCGG});
set(DBCCGG.nodeproperty_sel,'callback',{@DBCCGGnodeprop,DBCCGG});
set(DBCCGG.nodeproperty_listbox,'callback',{@DBCCGGnodeshow,DBCCGG});
set(DBCCGG.Degthr_hub,'callback',{@DBCCGGDeghub,DBCCGG});
set(DBCCGG.Degthr_othertxt,'callback',{@DBCCGGDegother,DBCCGG});
set(DBCCGG.Degthr_otheredi,'callback',{@DBCCGGDegotheredit,DBCCGG});
set(DBCCGG.Opt_show,'callback',{@DBCCGGshow,DBCCGG});
set(DBCCGG.Opt_exit,'callback',{@DBCCGGexit,DBCCGG});
end
function DBCCGGinputmat(varargin)
DBCCGG = varargin{3};
[fil pat ext] = uigetfile('*.txt','Connect Matrix');
set(DBCCGG.IO_inputedi,'string',fullfile(pat,fil));
end
function DBCCGGinputmod(varargin)
DBCCGG = varargin{3};
[fil pat ext] = uigetfile('*.txt','Select Module definition');
set(DBCCGG.Mod_inputedi,'string',fullfile(pat,fil));
Mod = load(fullfile(pat,fil));
DBCCGG.Mod = Mod;
MODNUM = length(unique(Mod));
set(DBCCGG.modc_txtnum,'string',num2str(MODNUM));
Modc = 1:MODNUM;
[pat,nam,ext] = fileparts(which('DynamicBC_CONGRAMGUI.m'));
load(fullfile(pat,'COLORMODULE.mat'));
COLORMODULEuse = COLORMODULE(1:MODNUM,:);
image(Modc,'parent',DBCCGG.modax,'CDataMapping','scaled');
colormap(DBCCGG.modax,COLORMODULEuse);axis(DBCCGG.modax,'off');
if isempty(dir(fullfile(pat,'COLTEMP.mat')))
    save(fullfile(pat,'COLTEMP.mat'),'COLORMODULEuse');
else
    delete(fullfile(pat,'COLTEMP.mat'))
    save(fullfile(pat,'COLTEMP.mat'),'COLORMODULEuse');
end    
end
function DBCCGGmodcol(varargin)
DBCCGG = varargin{3};
numCOL = get(DBCCGG.modc_txtnum,'string');
NUMCOL = str2num(numCOL);
[pat,nam,ext] = fileparts(which('DynamicBC_CONGRAMGUI.m'));
load(fullfile(pat,'COLORMODULE.mat'));
if isempty(dir(fullfile(pat,'COLTEMP.mat')))
    COLORMODULEuse = COLORMODULE(1:NUMCOL,:);
else
    load(fullfile(pat,'COLTEMP.mat'));    
end
selc.H = figure('units','norm','pos',[0.1 0.1 0.8 0.8],'name','Select Color');
indi = 1;
for j = 1:3
    for i = 1:6
        selc.ubg(indi) = uibuttongroup('parent',selc.H,'units','norm','pos',[(i-1)/6,0.9-j*0.8/3,1/6,1/3*0.8]);
        indi = indi+1;
    end
end
for i = 1:indi-1
    selc.ax(i) = axes('parent',selc.ubg(i),'units','norm','pos',[0,0.3,1,0.7]);
    image(i,'parent',selc.ax(i),'CDataMapping','scaled');colormap(selc.ax(i),COLORMODULE);
    set(selc.ax(i),'Clim',[1 18]);axis(selc.ax(i),'off');
%     selc.coledi(i) = uicontrol('parent',selc.ubg(i),'units','norm','pos',[0 0.1 1 0.1],'style','edit','string',num2str(COLORMODULE(i,:)));
    selc.coledi(i) = uicontrol('parent',selc.ubg(i),'units','norm','pos',[0 0.1 1 0.1],'style','edit','string',num2str(floor(COLORMODULE(i,:)*255))); % change to 0-255
    set(selc.coledi(i),'enable','off');
end
COLORMODULE_new = COLORMODULE;
COLORMODULE_new(1:NUMCOL,:) = COLORMODULEuse;
for i = 1:NUMCOL 
%     selc.ax(i) = axes('parent',selc.ubg(i),'units','norm','pos',[0,0.3,1,0.7]);
    image(i,'parent',selc.ax(i),'CDataMapping','scaled');colormap(selc.ax(i),COLORMODULE_new);
    set(selc.ax(i),'Clim',[1 18]);axis(selc.ax(i),'off');
%     set(selc.coledi(i),'string',num2str(COLORMODULE_new(i,:)));
    set(selc.coledi(i),'string',num2str(floor(COLORMODULE_new(i,:)*255))); % set the range 0-255
    set(selc.coledi(i),'enable','on');
end
selc.ok = uicontrol('parent',selc.H,'units','norm','pos',[0.4 0.01 0.2 0.08],'style','pushbutton','string','ok','callback',{@selcok,selc,pat,NUMCOL,DBCCGG});
set(selc.coledi,'callback',{@selccolc,selc,NUMCOL,pat,COLORMODULE_new,COLORMODULE,COLORMODULEuse});
end
function DBCCGGnodeprop(varargin)
DBCCGG = varargin{3};
numofproans = inputdlg('num of property','num of property',1,{'2'});
numofprop = str2num(numofproans{1});
for i = 1:numofprop
    [fil pat ext] = uigetfile('*.txt',['prop: ',num2str(i)]);
    liststr{i,1} = fullfile(pat,fil);
    pops = load(liststr{i,1});
    queans = questdlg('>mean+std?','>mean+std?','yes','no','all one','yes');
    if strcmp(queans,'yes')
        listpstr{i,1} = '>mean+std';
    elseif strcmp(queans,'no')
        meanpops = mean(pops);
        listpstr{i,1} = ['Other:largerthan:',num2str(meanpops),',MAX:',num2str(max(pops)),',MIN:',num2str(min(pops))];
    else
        listpstr{i,1} = 'Allone';
    end
end
set(DBCCGG.nodeproperty_listbox,'string',liststr);
set(DBCCGG.nodeproperty_listbox2,'string',listpstr);
end
function DBCCGGnodeshow(varargin)
DBCCGG = varargin{3};
listbox1 = get(DBCCGG.nodeproperty_listbox,'string');
listbox2 = get(DBCCGG.nodeproperty_listbox2,'string');
val = get(DBCCGG.nodeproperty_listbox,'value');
set(DBCCGG.nodeproperty_listbox2,'value',val);
stringinfo = listbox2{val};
if strcmp(stringinfo(1),'>')
    set(DBCCGG.Degthr_hub,'value',1);
    set(DBCCGG.Degthr_othertxt,'val',0);
    set(DBCCGG.Degthr_otheredi,'enable','off')
elseif strcmp(stringinfo(1),'O')
    set(DBCCGG.Degthr_hub,'value',0);
    set(DBCCGG.Degthr_othertxt,'val',1);
    positionsn1 = find(stringinfo=='n');
    positionsn2 = find(stringinfo==',');
    thrv = stringinfo(positionsn1+2:positionsn2(1)-1);
    set(DBCCGG.Degthr_otheredi,'string',thrv);
    set(DBCCGG.Degthr_otheredi,'enable','on');
else
    set(DBCCGG.Degthr_hub,'value',0);
    set(DBCCGG.Degthr_othertxt,'val',1);
    set(DBCCGG.Degthr_otheredi,'string','NULL');
    set(DBCCGG.Degthr_otheredi,'enable','off');
end
end
function DBCCGGDeghub(varargin)
DBCCGG = varargin{3};
set(DBCCGG.Degthr_hub,'value',1);
set(DBCCGG.Degthr_othertxt,'val',0);
set(DBCCGG.Degthr_otheredi,'enable','off');
val = get(DBCCGG.nodeproperty_listbox2,'val');
listbox2 = get(DBCCGG.nodeproperty_listbox2,'string');
listbox2{val} = '>mean+std';
set(DBCCGG.nodeproperty_listbox2,'string',listbox2);
end
function DBCCGGDegother(varargin)
DBCCGG = varargin{3};
set(DBCCGG.Degthr_hub,'value',0);
set(DBCCGG.Degthr_othertxt,'val',1);
allquest = questdlg('all one?','all one?','yes','no','no');

val = get(DBCCGG.nodeproperty_listbox,'value');
listbox2 = get(DBCCGG.nodeproperty_listbox2,'string');
listbox1 = get(DBCCGG.nodeproperty_listbox,'string');
strings = listbox2{val};

if strcmp(allquest,'yes')
    listbox2{val} = 'Allone';
    set(DBCCGG.nodeproperty_listbox2,'string',listbox2);
    set(DBCCGG.Degthr_otheredi,'string','NULL');
    set(DBCCGG.Degthr_otheredi,'enable','off');
else
    set(DBCCGG.Degthr_otheredi,'enable','on');
    if strcmp(strings(1),'O')
        positionsn1 = find(strings=='n');
        positionsn2 = find(strings==',');
        thrv = stringinfo(positionsn1+2:positionsn2(1)-1);
        set(DBCCGG.Degthr_otheredi,'string',thrv);
    else
        pops=load(listbox1{val});
        meanpops = mean(pops);
        listbox2{val} = ['Other:largerthan:',num2str(meanpops),',MAX:',num2str(max(pops)),',MIN:',num2str(min(pops))];
        set(DBCCGG.nodeproperty_listbox2,'string',listbox2);
        set(DBCCGG.Degthr_otheredi,'string',num2str(meanpops));
    end
end
end
function DBCCGGDegotheredit(varargin)
DBCCGG = varargin{3};
strinfo = get(DBCCGG.Degthr_otheredi,'string');
val = get(DBCCGG.nodeproperty_listbox2,'val');
list2 = get(DBCCGG.nodeproperty_listbox2,'string');
stringinfo = list2{val};
positionsn1 = find(stringinfo=='n');
positionsn2 = find(stringinfo==',');
thrv = stringinfo(positionsn1+2:positionsn2(1)-1);
list2{val} = ['Other:largerthan:',strinfo,stringinfo(positionsn2(1):end)];
set(DBCCGG.nodeproperty_listbox2,'string',list2);
end
function DBCCGGshow(varargin)
DBCCGG = varargin{3};
matrix = load(get(DBCCGG.IO_inputedi,'string'));
nsize = size(matrix,1);
TRILMAT = tril(ones(nsize),-1);
TRIRMAT = triu(ones(nsize),1);
if isempty(matrix(find(TRILMAT)))||isempty(matrix(find(TRIRMAT)))
    matrix = matrix+matrix';
end
Mod = load(get(DBCCGG.Mod_inputedi,'string'));
list1 = get(DBCCGG.nodeproperty_listbox,'string');
list2 = get(DBCCGG.nodeproperty_listbox2,'string');
Outliers = ones(nsize,length(list1));
for i = 1:length(list1)
    list2t = list2{i};
    Outlierst = load(list1{i});
    if strcmp(list2t(1),'A')
        Props(i) = NaN;
    elseif strcmp(list2t(1),'>')
        Outliers(:,i) = Outlierst;
        Props(i) = mean(Outlierst)+std(Outlierst);
    else
        Outliers(:,i) = Outlierst;
        positionsn1 = find(list2t=='n');
        positionsn2 = find(list2t==',');
        thrv = list2t(positionsn1+2:positionsn2(1)-1);
        Props(i) = str2num(thrv);
    end
end
numofcol = get(DBCCGG.modc_txtnum,'string');
NUMCOL = str2num(numofcol);

[pat,nam,ext] = fileparts(which('DynamicBC_CONGRAMGUI.m'));
load(fullfile(pat,'COLTEMP.mat'));
coloredge_2 = [0.8 0.8 0.8;1 0 0];
edgthr = get(DBCCGG.Edgethr_edi,'string');
DEGTHR = str2num(edgthr);
thrtype = get(DBCCGG.Edgethr_abs,'val');
if thrtype
    matrix = abs(matrix);
end
matrix(1:size(matrix,1)+1:end) = 0;
TCtype = get(DBCCGG.Edgethr_cost,'val');
if TCtype
    T = costthr(matrix,DEGTHR);
else
    T = DEGTHR;
end
Rmapthred = matrix.*(matrix>T);
Deg = Outliers(:,1);
[OutMatrix,OutMod,roire,ZeroInd,ReGroupInd,Outsortinfo] = DynamicBC_RegroupDat(Rmapthred,Mod,Deg);
rextend(:,1) = Outsortinfo;
for i0 = 2:length(list1)
    OutlierResortT = zeros(size(Outsortinfo));
    for i = 1:length(roire)
        if roire(i)~=0
            OutlierResortT(i) = Outliers(roire(i),i0);
        end
    end
    rextend(:,i0) = OutlierResortT;
end
for i = 1:length(OutMod)
    if OutMod(i)~=0
        ColorMod{i} = COLORMODULEuse(OutMod(i),:);
    else
        ColorMod{i} = [NaN NaN NaN];
    end
end
colorFrame = ColorMod;
for ip = 1:length(Props)
    if isnan(Props(ip))
        ColorMod1 = ColorMod;
    else
        DEGTHR = Props(ip);
        Outsortinfo = rextend(:,ip);
        for i = 1:length(OutMod)
            if OutMod(i)~=0
                if Outsortinfo(i)>DEGTHR
                    ColorMod1{i} = COLORMODULEuse(OutMod(i),:);
                else
                    ColorMod1{i} = [0.8 0.8 0.8];
                end
            else
                ColorMod1{i} = [NaN NaN NaN];
            end
        end
    end
    colorFrame2{ip} = ColorMod1;
end
linewid = get(DBCCGG.Edgewidth_equ,'val');
if linewid
    line_widthtype = 1;
else 
    line_widthtype = 2;
end
linfac = get(DBCCGG.Edgewidth_chaedi,'string');
linefactor = str2num(linfac);
abov = get(DBCCGG.Overlap_with,'val');
if abov
    abovetype = 1;
else
    abovetype = 2;
end
colt = get(DBCCGG.Edgeshowtype_gray,'val');
if colt
    colortype = 1;
else
    colortype = 2;
end
modcolor = COLORMODULEuse;
nsizet = get(DBCCGG.Nodsize_edi,'string');
nodes = str2num(nsizet);
numSlices = length(OutMod);
% save temps
DynamicBC_connectogramOutlier(numSlices,rextend,colorFrame2,colorFrame,OutMod,OutMatrix,line_widthtype,abovetype,linefactor,colortype,modcolor,nodes)
DynamicBC_connectogramOutlierA(numSlices,rextend,colorFrame2,colorFrame,OutMod,OutMatrix,line_widthtype,abovetype,linefactor,colortype,modcolor,nodes)
DynamicBC_connectogramOutlierB(numSlices,rextend,colorFrame2,colorFrame,OutMod,OutMatrix,line_widthtype,abovetype,linefactor,colortype,modcolor,nodes)
DynamicBC_connectogramOutlierW(numSlices,rextend,colorFrame2,colorFrame,OutMod,OutMatrix,line_widthtype,abovetype,linefactor,colortype,modcolor,nodes)
end
function DBCCGGexit(varargin)
DBCCGG = varargin{3};
close(DBCCGG.fig);
end
function selccolc(varargin)
selc = varargin{3};
NUMCOL = varargin{4};
pat = varargin{5};
COLORMODULE_new = varargin{6};
COLORMODULE = varargin{7};
COLORMODULEuse = varargin{8};
for i = 1:NUMCOL
    COLS = get(selc.coledi(i),'string');
    COLORS = str2num(COLS);
%     COLORMODULE_new(i,:) = COLORS;
    COLORMODULE_new(i,:) = COLORS/255; % set color range 0-255
end
COLORMODULEuse(1:NUMCOL,:) = COLORMODULE_new(1:NUMCOL,:);
save(fullfile(pat,'COLTEMP.mat'),'COLORMODULEuse');
for i = 1:NUMCOL
    image(i,'parent',selc.ax(i),'CDataMapping','scaled');colormap(selc.ax(i),COLORMODULE_new);
    set(selc.ax(i),'Clim',[1 18]);axis(selc.ax(i),'off');
%     set(selc.coledi(i),'string',num2str(COLORMODULE_new(i,:)));
%     set(selc.coledi(i),'enable','on');
end
end
function selcok(varargin)
selc = varargin{3};
pat = varargin{4};
NUMCOL = varargin{5};
DBCCGG = varargin{6};
close(selc.H)
load(fullfile(pat,'COLTEMP.mat'));
Modc = 1:NUMCOL;
image(Modc,'parent',DBCCGG.modax,'CDataMapping','scaled');colormap(DBCCGG.modax,COLORMODULEuse);axis(DBCCGG.modax,'off');
end
function T=costthr(matrix,edgthr)
defcost = 0;
steps = 0.0001;
T=max(matrix(:));
while defcost<edgthr
    mat = matrix>T;
    defcost = nnz(mat)/(size(mat,1)*(size(mat,1)-1));
    T = T-steps;
end
T = T+steps;
end
