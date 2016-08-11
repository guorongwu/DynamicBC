function [M] = DynamicBC_PWGC_SLW(data,nvar,flag_1n, flag_gcd, pvalue,order,seed_signal,flag_par,m)
% Compute correlation matrix A.
warning off
disp('abc')
if nargin<8
    m = 5000; %block size
    flag_par = 0; % open parllel workers
elseif nargin==8
    m = 5000;
end
nobs = size(data,1);
M.nvar=nvar;
M.nobs=nobs;
M.pvalue=pvalue;
gc_th = wgr_inv_pwGC_F(pvalue,nobs,order);
M.gc_th=gc_th;

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
        RR = wgr_pwGC_OLS(X,seed_signal,order,flag_par);
        M.matrix{i,1} = RR;
    else  % FC: n x n
        for j=1:nbin
            if j~=nbin
                ind_Y = (j-1)*m+1:j*m ;
            else
                ind_Y = (j-1)*m+1:nvar ;
            end

            Y=data(:,ind_Y);            
            RR = wgr_pwGC_OLS(X,Y,order,flag_par);
            if flag_gcd
                M.matrix{i,j} = RR.*(RR>=gc_th);
            else
                M.matrix{i,j} = RR;
            end
        end
    end
end  


function [pwgc] = wgr_pwGC_OLS(data1,data2,order,flag_par);

[nvar1]=size(data1,2);
[nvar2]=size(data2,2);

if ~nvar2
    pwgc=zeros(nvar1,nvar1);
    if flag_par
        parfor drive=1:nvar1
            for target=1:nvar1
                if drive ~= target
                    pwgc(drive,target) = wgr_GCA_OLS(data1(:,drive),data1(:,target),[],order); 
                end
            end
        end
    else
        for drive=1:nvar1
            for target=1:nvar1
                if drive ~= target
                    pwgc(drive,target) = wgr_GCA_OLS(data1(:,drive),data1(:,target),[],order); 
                end
            end
        end
    end
else 
    if nvar2==1;
        pwgc=zeros(2,nvar1);
        if flag_par
            parfor target=1:nvar1
                pwgc(1,target) = wgr_GCA_OLS(data1(:,target),data2,[],order); 
            end
            parfor target=1:nvar1
                pwgc(2,target) = wgr_GCA_OLS(data2,data1(:,target),[],order); 
            end
        else
            for target=1:nvar1
                pwgc(1,target) = wgr_GCA_OLS(data1(:,target),data2,[],order); 
                pwgc(2,target) = wgr_GCA_OLS(data2,data1(:,target),[],order); 
            end
        end    
    else
        pwgc=zeros(nvar1,nvar2);
        if flag_par
            parfor drive=1:nvar1
                for target=1:nvar2
                    pwgc(drive,target) = wgr_GCA_OLS(data1(:,drive),data2(:,target),[],order); 
                end
            end
        else
            for drive=1:nvar1
                for target=1:nvar2
                    pwgc(drive,target) = wgr_GCA_OLS(data1(:,drive),data2(:,target),[],order); 
                end
            end
        end
    end        
end

function gc = wgr_inv_pwGC_F(pvalue,nobs,porder);
n = nobs - porder;
th = 1+finv(1-pvalue,porder,n-2*porder-1)/(n-2*porder-1)*porder;
gc = log(th);

function [cgc] = wgr_GCA_OLS(y,x,z,order)
%% y->x condtion on z. based on covariance matrix
%% x: N*nx, y: N*ny, z:N*nz.
%% 
[N,nx]=size(x);
%now
X = x(order+1:end,:);
%past
past_ind = repmat([1:order],N-order,1) + repmat([0:N-order-1]',1,order);
xz = [x z];
XZ_past = reshape(xz(past_ind,:),N-order,order*size(xz,2));
Y_past = reshape(y(past_ind,:),N-order,order*size(y,2));
XZY_past = [XZ_past Y_past];

% Remove mean
xzyc = bsxfun(@minus,XZY_past,sum(XZY_past,1)/(N-order)); 
Xc = bsxfun(@minus,X,sum(X,1)/(N-order)); 
xzc = bsxfun(@minus,XZ_past,sum(XZ_past,1)/(N-order)); 

%Covariance matrix
cov_X = (Xc' * Xc) / (N-order-1); 
cov_xz = (xzc' * xzc) / (N-order-1);
cov_xzy = (xzyc' * xzyc) / (N-order-1);
cov_X_xz = (Xc' * xzc) / (N-order-1); 
cov_X_xzy = (Xc' * xzyc) / (N-order-1);

%Partial cross-covariance
cov_X_xz = cov_X - cov_X_xz/cov_xz*cov_X_xz';
cov_X_xyz = cov_X - cov_X_xzy/cov_xzy*cov_X_xzy';
cgc = log(det(cov_X_xz)/det(cov_X_xyz));