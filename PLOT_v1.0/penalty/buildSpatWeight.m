function R = buildSpatWeight(R, w, G, Gopt, ni, yi)
%
% Build the spatially variant weighting factor for regularizer
%
% Guobao Wang @ UC Davis (10-01-2012)
%

% check
siz = R.imgsiz;
if nargin<2 | isempty(w)
    w = ones(prod(siz),1);
end

% caring sensitivity
if nargin>2
    w1 = proj_back(G, Gopt, ni);
    wy = proj_back(G, Gopt, ni.*yi);
    wx = w1./max(mean(wy)*1e-9,wy);
    wx(w1==0&wy==0) = 1;
    w  = w.*wx;
end

% set the central pixel to 1
w = w / w(sub2ind(siz,round(siz(1)/2),round(siz(2)/2),round(siz(3)/2)));

% spatial weight
R.spatwgt = w(:);

