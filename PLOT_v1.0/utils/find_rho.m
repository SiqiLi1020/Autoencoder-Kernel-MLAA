function [rho_opt, kap] = find_rho(wi, rho)
%
% use circulant approximation to find an empirically good parameter rho for 
% diagonal matrix A = D[w] + rho*I that appears in ADMM algorithms
%
% gbwang@ucdavis.edu (01-29-2013)
%

for i = 1:length(rho)
    kap(i) = max(wi+rho(i))./min(wi+rho(i));
end

% choose
[a,i] = max(abs(gradient(kap)));
rho_opt = rho(i);
