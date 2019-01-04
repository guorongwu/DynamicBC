function DynamicBC_clustermaps(num_clusters,outputd,Maskimg,subjdir,dmethod,filter_str,flag_estimate)
outputd = strcat(outputd);
subjdir = strcat(subjdir);
SubFold = dir(subjdir); SubFold(1:2)=[];
if isempty(SubFold)
    error(sprintf('No subject data in %s\n',subjdir))
    return;
else
    tmp=[];
    for i=1:length(SubFold) %check for mac os
        if ~SubFold(i).isdir
            tmp = [tmp i];
        end
    end
    SubFold(tmp)=[]; clear tmp
end
NumOfSubFold = length(SubFold);
filter_str2 = strrep(strrep(filter_str,'*',''),'_','');
if isempty(filter_str2)
    filter_str3='';
else
    filter_str3 = ['_',filter_str2];
end
dirtemp = fullfile(subjdir,SubFold(1).name);
[v] = wgr_read_dir_NIFTI(dirtemp,filter_str);
vm = spm_vol(Maskimg);
if any(vm(1).dim-v(1).dim)
    fprintf('Data DIM:[%d,%d,%d]\n',v(1).dim)
    fprintf('Mask %s, DIM:[%d,%d,%d]\n',Maskimg,vm(1).dim)
    fprintf('(Nearest Neighbour)Reslice mask ... \n')
    [fpath_seed,name,ext] = fileparts(Maskimg);
    tmpext = strfind(ext,',');
    if ~isempty(tmpext)
        ext = ext(1:tmpext(1)-1);
    end
    Maskimg0 = fullfile(fpath_seed,[name,'_resliced',ext]);
    dynamicBC_Reslice(Maskimg,Maskimg0,0,v(1).fname);
    vm = spm_vol(Maskimg0);
end
MASKdat = spm_read_vols(vm);
MASKdat(isnan(MASKdat))=0;
maskind = find(MASKdat);
dims = vm(1).dim;

alen = 0; Subindex={}; 
DATA={};  DATAmx={};
for i = 1:NumOfSubFold
    dirtemp = fullfile(subjdir,SubFold(i).name);
    [~,Dat] = wgr_read_dir_NIFTI(dirtemp,filter_str);
    Dat(isnan(Dat)) = 0; %nan
    blen = alen+size(Dat,2);
    Subindex{i,1} = alen+1:blen;
    tmp = Dat(maskind,:);
    DATA{1,i} = tmp;
    tmp2 = std(tmp);
    id_max = DynamicBC_extrema(tmp2);
    Subindex{i,2} = id_max;
    DATAmx{1,i} = tmp(:, id_max);
    alen = blen;
end

DATAmx = cell2mat(DATAmx);
DATA = cell2mat(DATA);
%%
num_clusters_max = max(num_clusters);
Repeats= 5; 

num_clustersall = [];
if flag_estimate
    fprintf('The maximum number of data clusters = %d\n',num_clusters_max)
    IDXest = zeros(size(DATAmx,2),num_clusters_max);
    distortion= [];
    for i=1:num_clusters_max
        [IDXest(:,i),Cest{i},sumdest{i}] = kmeans(DATAmx',i,'emptyaction','singleton','replicate',Repeats, 'empty', 'drop');
        distortion(i,1) = sum(sumdest{i});  
    end
    variance=distortion(1:end-1)-distortion(2:end);
    distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
%     [r,~]=find(distortion_percent>0.9);
%     K_elbow = r(1,1)+1;
    
    crit = {'silhouette','CalinskiHarabasz', 'DaviesBouldin'};    %     crit = {'CalinskiHarabasz', 'DaviesBouldin', 'gap', 'silhouette'};
    fig = figure('color','w','units','norm','pos',[0.1,0.3,0.7,0.25]);
    
    %https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set
    subplot(1,4,1); plot(2:num_clusters_max,distortion_percent,'b-o')
    xlabel('Number of Clusters'); 
    ylabel('Percent of variance explained') %Percentage of variance explained is the ratio of the between-group variance to the total variance,
    OptimalK=[];
    for i=1:length(crit)
        eva{i} = evalclusters(DATAmx',IDXest,crit{i});
        subplot(1,4,i+1); plot(eva{i})%       
        tmp=[];
        tmp.OptimalK = eva{i}.OptimalK;
        OptimalK(i) = eva{i}.OptimalK;
        fprintf('(%s) Optimal K=%d\n',crit{i},OptimalK(i))
        tmp.OptimalY = eva{i}.OptimalY;
        tmp.InspectedK = eva{i}.InspectedK;
        tmp.CriterionValues = eva{i}.CriterionValues; 
        tmp.CriterionName = eva{i}.CriterionName;
        tmp.ClusteringFunction = eva{i}.ClusteringFunction;
        tmp.NumObservations = eva{i}.NumObservations;
        tmp.Missing = eva{i}.Missing;        
        eva{i} = tmp;
    end
    
    tmp  = ceil(mean(OptimalK));
    fprintf('(Mean) Optimal K=%d\n',tmp)
    num_clustersall = [num_clustersall tmp];
    if length(num_clusters)>1
        tmp = num_clusters;  tmp  = sort(tmp); tmp(end)=[];
        num_clustersall = [num_clustersall tmp];
    end
    if ~isempty(filter_str2)
        set(fig,'name',filter_str2)
        saveas(fig,fullfile(outputd,['EVA','_',filter_str2,'.fig']))
        save(fullfile(outputd,['EVA','_',filter_str2,'.mat']),'eva','num_clusters','OptimalK')
    else        
        saveas(fig,fullfile(outputd,'EVA.fig'))
        save(fullfile(outputd,'EVA.mat'),'eva','num_clusters','OptimalK')
    end
else
    num_clustersall = [num_clusters];
end

num_clustersall  = unique(num_clustersall);
for inum_clusters = num_clustersall
    fprintf('K-means Clustering: K=%d\n',inum_clusters)
    if ~flag_estimate
        [IDX,C,sumd] = kmeans(DATAmx', inum_clusters, 'distance', dmethod, 'Replicates', Repeats, 'empty', 'drop');
    else
        C = Cest{inum_clusters};
    end
    [IDXall, Call, SUMDall, Dall] = kmeans(DATA', inum_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'empty', 'drop', 'Start', C);
    if ~isempty(filter_str2)
        save_dir = fullfile(outputd,[dmethod,'_Kmeans_',num2str(inum_clusters),'_',filter_str2]);
    else
        save_dir = fullfile(outputd,[dmethod,'_Kmeans_',num2str(inum_clusters)]);
    end
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    if inum_clusters==num_clustersall(1)
        v1 = vm(1);
        v1.dt = [16,0];
        v1.n=[1,1];
    end

    for i = 1:inum_clusters
        DAT = zeros(dims);
        DAT(maskind) = C(i,:);
        if ~isempty(filter_str2)
            v1.fname = fullfile(save_dir,['Centroid_',num2str(i),'_',filter_str2,'.nii']); 
        else
            v1.fname = fullfile(save_dir,['Centroid_',num2str(i),'.nii']); 
        end
        spm_write_vol(v1, DAT);
    end

    %subject's median maps for each cluster
    IDX_subj={};
    for isubj = 1:NumOfSubFold
        fprintf('.')
        for i = 1:inum_clusters
            indtemp = find(IDXall(Subindex{isubj,1})==i);
            if ~isempty(indtemp)
                DAT(maskind) = nanmedian(DATA(:,indtemp),2);
                subdir = fullfile(save_dir,SubFold(isubj).name); 
                if ~exist(subdir,'dir')
                    mkdir(subdir);
                end
                if ~isempty(filter_str2)
                    v1.fname = fullfile(subdir,['MedianCluster',num2str(i),'_',filter_str2,'.nii']);
                else
                    v1.fname = fullfile(subdir,['MedianCluster',num2str(i),'.nii']);
                end
                spm_write_vol(v1,DAT);
            end
        end
        IDX_subj{isubj} = IDXall(Subindex{isubj,1});
    end
    K = inum_clusters;
    save(fullfile(save_dir,['cluster_index',filter_str3,'.mat']),'IDX_subj','K','Subindex','IDXall','SubFold')
    fprintf('\n')
end
clear DATAmx DATA
end

function [v,data] = wgr_read_dir_NIFTI(sub_dir,filter_str)
niiname = dir(fullfile(sub_dir,[filter_str,'*.nii']));
if isempty(niiname)
    niiname = dir(fullfile(sub_dir,[filter_str,'*.img']));
    if isempty(niiname)             
        warning(['Warning: no NIFTI file exist in: ',sub_dir]);
    end
end

num_volume = length(niiname) ;
if num_volume==1
    disp(['Only one file in: ',sub_dir])
    niifile = fullfile(sub_dir,niiname.name);
    v = spm_vol(niifile);
    data = spm_read_vols(v);
    nobs = size(data,4);
    data = reshape(data,[],nobs); %nvar*nobs
else
    nobs = num_volume;
    tmp_dat = spm_read_vols(spm_vol(fullfile(sub_dir,niiname(1).name)));
    data = zeros(numel(tmp_dat),nobs); %nvar*nobs
    for j=1:nobs
        niifile = fullfile(sub_dir,niiname(j).name);
        v = spm_vol(niifile);
        tmp_dat = spm_read_vols(v);
        data(:,j) = tmp_dat(:);
    end
end

end

%% subfunction
function dynamicBC_Reslice(PI,PO,hld,TargetSpace)
%   PI - input filename
%   PO - output filename
%   NewVoxSize - 1x3 matrix of new vox size.
%   hld - interpolation method. 0: Nearest Neighbour. 1: Trilinear.
%   TargetSpace - Define the target space. 'ImageItself': defined by the image itself (corresponds  to the new voxel size); 'XXX.img': defined by a target image 'XXX.img' (the NewVoxSize parameter will be discarded in such a case).

headIN = spm_vol(TargetSpace) ;
headIN = headIN(1);
dataIN = spm_read_vols(headIN);
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
end
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
end

function [imax] = DynamicBC_extrema(x)
%EXTREMA   Gets the global extrema points from a time series.
%   [XMAX,IMAX,XMIN,IMIN] = EXTREMA(X) returns the global minima and maxima 
%   points of the vector X ignoring NaN's, where
%    XMAX - maxima points in descending order
%    IMAX - indexes of the XMAX
%    XMIN - minima points in descending order
%    IMIN - indexes of the XMIN
%
%   DEFINITION (from http://en.wikipedia.org/wiki/Maxima_and_minima):
%   In mathematics, maxima and minima, also known as extrema, are points in
%   the domain of a function at which the function takes a largest value
%   (maximum) or smallest value (minimum), either within a given
%   neighbourhood (local extrema) or on the function domain in its entirety
%   (global extrema).
%
%   Example:
%      x = 2*pi*linspace(-1,1);
%      y = cos(x) - 0.5 + 0.5*rand(size(x)); y(40:45) = 1.85; y(50:53)=NaN;
%      [ymax,imax,ymin,imin] = extrema(y);
%      plot(x,y,x(imax),ymax,'g.',x(imin),ymin,'r.')
%
%   See also EXTREMA2, MAX, MIN

%   Written by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico, 2004
%
%   nubeobscura@hotmail.com

% From       : http://www.mathworks.com/matlabcentral/fileexchange
% File ID    : 12275
% Submited at: 2006-09-14
% 2006-11-11 : English translation from spanish. 
% 2006-11-17 : Accept NaN's.
% 2007-04-09 : Change name to MAXIMA, and definition added.


xmax = [];
imax = [];
xmin = [];
imin = [];

% Vector input?
Nt = numel(x);
if Nt ~= length(x)
 error('Entry must be a vector.')
end

% NaN's:
inan = find(isnan(x));
indx = 1:Nt;
if ~isempty(inan)
 indx(inan) = [];
 x(inan) = [];
 Nt = length(x);
end

% Difference between subsequent elements:
dx = diff(x);

% Is an horizontal line?
if ~any(dx)
 return
end

% Flat peaks? Put the middle element:
a = find(dx~=0);              % Indexes where x changes
lm = find(diff(a)~=1) + 1;    % Indexes where a do not changes
d = a(lm) - a(lm-1);          % Number of elements in the flat peak
a(lm) = a(lm) - floor(d/2);   % Save middle elements
a(end+1) = Nt;

% Peaks?
xa  = x(a);             % Serie without flat peaks
b = (diff(xa) > 0);     % 1  =>  positive slopes (minima begin)  
                        % 0  =>  negative slopes (maxima begin)
xb  = diff(b);          % -1 =>  maxima indexes (but one) 
                        % +1 =>  minima indexes (but one)
imax = find(xb == -1) + 1; % maxima indexes
imin = find(xb == +1) + 1; % minima indexes
imax = a(imax);
imin = a(imin);

nmaxi = length(imax);
nmini = length(imin);                

% Maximum or minumim on a flat peak at the ends?
if (nmaxi==0) && (nmini==0)
 if x(1) > x(Nt)
  xmax = x(1);
  imax = indx(1);
  xmin = x(Nt);
  imin = indx(Nt);
 elseif x(1) < x(Nt)
  xmax = x(Nt);
  imax = indx(Nt);
  xmin = x(1);
  imin = indx(1);
 end
 return
end

% Maximum or minumim at the ends?
if (nmaxi==0) 
 imax(1:2) = [1 Nt];
elseif (nmini==0)
 imin(1:2) = [1 Nt];
else
 if imax(1) < imin(1)
  imin(2:nmini+1) = imin;
  imin(1) = 1;
 else
  imax(2:nmaxi+1) = imax;
  imax(1) = 1;
 end
 if imax(end) > imin(end)
  imin(end+1) = Nt;
 else
  imax(end+1) = Nt;
 end
end
xmax = x(imax);
% xmin = x(imin);

% NaN's:
if ~isempty(inan)
 imax = indx(imax);
%  imin = indx(imin);
end

% Same size as x:
imax = reshape(imax,size(xmax));
% imin = reshape(imin,size(xmin));

% % Descending order:
% [temp,inmax] = sort(-xmax); clear temp
% xmax = xmax(inmax);
% imax = imax(inmax);
% [xmin,inmin] = sort(xmin);
% imin = imin(inmin);


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com
end