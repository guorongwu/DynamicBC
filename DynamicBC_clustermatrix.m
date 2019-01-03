function [save_dir,IDX_subj]=DynamicBC_clustermatrix(num_clusters,outputd,matname,subjdir,dmethod,flag_estimate)
evalstr = ['mats = ',matname,';'];
SubFold = dir(subjdir); SubFold(1:2)=[];
if isempty(SubFold)
    error(sprintf('No subject data in %s\n',subjdir))
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
NumOfSubFold = size(SubFold,1);
alen = 0; 
inums = 1; Subindex={}; 
for i = 1:NumOfSubFold
    sub_dir = fullfile(subjdir,SubFold(i).name);
    if isdir(sub_dir)
        matsnames = dir(fullfile(sub_dir,'*.mat'));
        load(fullfile(sub_dir,matsnames.name))
        eval(evalstr);
        if i==1
            tt = mats{1}-mats{1}'; tt(tt<10^-6)=0;
            if any(tt(:))
                id = find([ones(size(tt))-eye(size(tt))]);
                disp('asymmetry matrix')
            else
                id = find(tril(ones(size(tt)),-1));
                disp('symmetry matrix')
            end
            clear tt  
        end
        nmats = length(mats);
        blen = alen+nmats;
        Subindex{i,1} = alen+1:blen;
        for imat = 1:nmats
            matrixs(:,imat) = reshape(mats{imat},numel(mats{imat}),1);
        end
        matrixs = full(matrixs);
        matrix_all{1,inums} = matrixs;
        subjlist{inums} = SubFold(i).name;
        tmp2 = std(matrixs);
        id_max = DynamicBC_extrema(tmp2);
        Subindex{i,2} = id_max;
        DATAmx{1,i} = matrixs(:, id_max);
        alen = blen;
        inums = inums+1;
    end
end
dims = size(mats{1});
matrixp = cell2mat(DATAmx);
matrixall =  cell2mat(matrix_all);
%%
num_clusters_max = max(num_clusters);
Repeats=10; 

num_clustersall = [];
if flag_estimate
    fprintf('The maximum number of data clusters = %d\n',num_clusters_max)
    IDXest = zeros(size(matrixp,2),num_clusters_max);
    distortion= [];
    for i=1:num_clusters_max
        [IDXest(:,i),Cest{i},sumdest{i}] = kmeans(matrixp',i,'emptyaction','singleton','replicate',Repeats, 'empty', 'drop');
        distortion(i,1) = sum(sumdest{i});  
    end
    variance=distortion(1:end-1)-distortion(2:end);
    distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
    
    crit = {'silhouette','CalinskiHarabasz', 'DaviesBouldin'};    %     crit = {'CalinskiHarabasz', 'DaviesBouldin', 'gap', 'silhouette'};
    fig = figure('color','w','units','norm','pos',[0.3,0.3,0.6,0.3]);
    %https://en.wikipedia.org/wiki/Determining_the_number_of_clusters_in_a_data_set
    subplot(1,4,1); plot(2:num_clusters_max,distortion_percent,'b-o')
    xlabel('Number of Clusters'); 
    ylabel('Percent of variance explained') %Percentage of variance explained is the ratio of the between-group variance to the total variance,
    OptimalK=[];
    for i=1:length(crit)
        eva{i} = evalclusters(matrixp',IDXest,crit{i});
        subplot(1,4,i+1); plot(eva{i})%         subplot(2,2,i); plot(eva{i})
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
    saveas(fig,fullfile(outputd,'EVA.fig'))
    save(fullfile(outputd,'EVA.mat'),'eva','num_clusters','OptimalK')
else
    num_clustersall = [num_clusters];
end

num_clustersall  = unique(num_clustersall);
for inum_clusters = num_clustersall
    fprintf('K-means Clustering: K=%d\n',inum_clusters)
    if ~flag_estimate
        [IDX,C,sumd,D] = kmeans(matrixp', inum_clusters, 'distance', dmethod, 'Replicates', Repeats, 'empty', 'drop');
    else
        C = Cest{inum_clusters};
    end
    [IDXall, Call, SUMDall, Dall] = kmeans(matrixall', inum_clusters, 'distance', dmethod, 'Replicates', 1, 'Display', 'iter', 'empty', 'drop', 'Start', C);
    save_dir = fullfile(outputd,['mat_',dmethod,'_Kmeans_',num2str(inum_clusters)]);
    if ~exist(save_dir,'dir')
        mkdir(save_dir);
    end
    for i = 1:inum_clusters
        DAT = reshape(C(i,:),dims(1),dims(2));
        save(fullfile(save_dir,['pCentroid_',num2str(i),'.mat']),'DAT')
        DAT = reshape(Call(i,:),dims(1),dims(2));
        save(fullfile(save_dir,['Centroid_',num2str(i),'.mat']),'DAT')
    end

    IDX_subj={};
    for isubj = 1:inums-1
        matrix = matrix_all{isubj};
        subdir = fullfile(save_dir,subjlist{isubj});
        if ~exist(subdir,'dir')
            mkdir(subdir);
        end
        for i = 1:inum_clusters
            indtemp = find(IDXall(Subindex{isubj,1})==i);
            meanmaptemp = nanmedian(matrix(:,indtemp),2);
            DAT = reshape(meanmaptemp,dims(1),dims(2));
            save(fullfile(subdir,['MedianCluster',num2str(i),'.mat']),'DAT')
        end
        IDX_subj{isubj} = IDXall(Subindex{isubj,1});
    end
    K = inum_clusters;
    save(fullfile(save_dir,['cluster_index.mat']),'IDX_subj','K','Subindex','IDXall','SubFold')
    fprintf('\n')
end
clear matrixp matrixall matrix_all

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