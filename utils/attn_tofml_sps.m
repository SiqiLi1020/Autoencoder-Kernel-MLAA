function [x, li, L] = attn_ml_sps(yt, bt, G, Gopt, x, rt, maxit, li, aa)
% This is a sub program for attenuation map reconstruction from PET emssion
% data
%
% gbwang@ucdavis.edu 05-01-2018
%

%% check inputs

numprj = prod(Gopt.prjsiz);
numfrm = size(yt,2);
numbin = size(yt,1)/numprj;

if nargin<2 | isempty(bt)
    bt = ones(numprj,numfrm);
end
if nargin<6 | isempty(rt)
    rt = ones(numprj*numbin,numfrm);
    yt = yt + rt;
end
if nargin<8
    li = [];
end
if nargin<9 
    aa = [];
end
if isempty(aa)
    aa = proj_forw(G, Gopt, ones(size(x)));
end


%% iterate
for it = 1:maxit
    
    % foward projection
    if isempty(li)
        li = proj_forw(G, Gopt, x);
    end
    
    % gradient
    lt = repmat(li,[1 numbin]);
    lt = repmat(lt(:),[1 numfrm]);
    yb = bt.*exp(-lt) + rt;
    yr = (1-yt./yb).*(yb-rt);
    yr = reshape(yr, [numprj numbin numfrm]);
    hi = sum(sum(yr,3),2);
    gx = proj_back(G, Gopt, hi);
    
    % optimal curvature and the optimization transfer weight
    nt = sum(trl_curvature(yt, bt, rt, lt, 'oc'),2);
    nt = reshape(nt,[numprj numbin]);
    wx = proj_back(G, Gopt, sum(nt,2).*aa, numfrm);
    
    % objective function
    L(it) = sum(yt(:).*log(yb(:)) - yb(:));
    
    % update image
    x(Gopt.mask) = x(Gopt.mask) + gx(Gopt.mask) ./ wx(Gopt.mask);
    x(Gopt.mask(:)&wx(:)==0) = 0;
    x = max(0,x);
    
    % update projection
    if it<maxit
        li = proj_forw(G, Gopt, x, numfrm);
    end
    
end
