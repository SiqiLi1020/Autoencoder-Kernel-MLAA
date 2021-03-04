function [xt, li, L] = emis_ml_em(yi, ci, G, Gopt, xt, ri, maxit, li, gg)
% This is a sub program for emission image reconstruction from PET data
%
% gbwang@ucdavis.edu 05-01-2018
%


%% check inputs
numprj = size(yi,1);
numfrm = size(yi,2);
numpix = prod(Gopt.imgsiz);
if nargin<2 | isempty(ci)
    ci = ones(numprj,1);
end
if size(ci,2)<numfrm
    ci = repmat(ci,[1 numfrm]);
end
if nargin<5 | isempty(xt)
    xt = ones(numpix,numfrm);
end
xt(xt<=0) = 1e-4*mean(xt(:));
if nargin<6 | isempty(ri)
    ri = zeros(numprj,numfrm);
end
if nargin<8 | isempty(li)
    li = proj_forw(G, Gopt, xt, numfrm);
end
if nargin<9 | isempty(gg)
    gg = proj_back(G, Gopt, ci, numfrm);
end
if size(gg,2)<numfrm
    gg = repmat(gg,[1 numfrm]);
end

% mask
mask = gg>0;
xt(~mask) = 0;

%% iterate
for it = 1:maxit
    
    % gradient
    yb = ci.*li + ri;
    yr = yi./yb;
    xb = proj_back(G, Gopt, ci.*yr, numfrm);
    xt = xt./gg .* xb;
    xt(gg==0) = 0;
    
    % objective function
    L(it) = sum(yi(:).*log(yb(:)) - yb(:));
    
    % update projection
    if it<maxit
        li = proj_forw(G, Gopt, xt, numfrm);
    end
    
end

