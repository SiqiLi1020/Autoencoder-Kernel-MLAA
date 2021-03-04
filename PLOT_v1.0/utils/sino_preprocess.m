function [yi, ri, ni] = sino_preprocess(yi, ri, ni)
% 
% pre-processing sinograms
%
% gbwang@ucdavis.edu (01-09-2013)
%

if isempty(ri)
    ri = zeros(size(yi));
end
if isempty(ni)
    ni = ones(size(yi));
end

% shift
minri = min(ri);
if minri<0
    ri = ri - minri;
    yi = yi - minri;
end
yi(yi<0) = 0;


% vectors
yi = yi(:);
ri = ri(:);
ni = ni(:);
