function [nuv_opt, nuw_opt, kap] = find_nuv(wi, ni, G, Gopt, nuv, wk, C, nuw)
%
% use circulant approximation to find an empirically good parameter nuv for
% matrix A = G'D[wi]G + nuv*C'D[wk]C + nuw*I that appear in ADMM algorithms
%
% gbwang@ucdavis.edu (01-29-2013)
%
% 

% check
if isempty(wi) 
    if isempty(Gopt) | ~isfield(Gopt,'prjsiz')
        wi = ones(size(G,1),1);
    else
        wi = ones(prod(Gopt.prjsiz),1);
    end
end
if isempty(wk) & not(isempty(C))
    wk = ones(size(C,1),1);
end
if isempty(ni)
    ni = ones(size(wi));
end
if nargin<8 | isempty(nuw)
    nuw = 0;
end

% cone coefficients
for i = 1:length(nuv)
    for j = 1:length(nuw)
        lam = cone_coef(wi, G, Gopt, ni, nuv(i), wk, C, nuw(j));
        lam = lam(lam(:)>0);
        kap(i,j) = max(lam)/min(lam);
    end
end

% choose
[a,k] = min(kap(:));
[i,j] = ind2sub([length(nuv),length(nuw)], k);
nuv_opt = nuv(i);
nuw_opt = nuw(j);
