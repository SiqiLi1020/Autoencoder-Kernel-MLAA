function lam = cone_coef(wy, G, Gopt, ni, alpha, w, C, beta)
%--------------------------------------------------------------------------
% calculate the fourier coefficients of the cone filter
%
% gbwang@ucdavis.edu, 02-01-2013
%

% check
imgsiz = Gopt.imgsiz;
if length(imgsiz)==2
    imgsiz = [imgsiz 1];
end
if alpha>0 & isempty(w)
    w = ones(size(C,1),1);
end
if nargin<8 | isempty(beta)
    beta = 0;
end

% make an image with value at the center
J = [ceil(imgsiz(1)/2),ceil(imgsiz(2)/2),ceil(imgsiz(3)/2)];
j = sub2ind(imgsiz, J(1), J(2), J(3));
e = zeros(imgsiz);
e(j) = 1;

% the central column of the Fisher information matrix
y = ni.*proj_forw(G, Gopt, e(:));
x = proj_back(G, Gopt, ni.*wy.*y);
if alpha>0
    x = x + alpha*(C'*(w.*(C*e(:))));
end
if beta>0
    x = x + beta*e(:);
end
x = reshape(x,imgsiz);
x = ShiftSym(x, imgsiz, J);

% the FFT coefficients
x_fft = fftn(fftshift(x));
lam = max(0,real(x_fft));


%--------------------------------------------------------------------------
function imgout = ShiftSym(img, imgsiz, MaxValPos)
%--------------------------------------------------------------------------

img_tmp = zeros(2*ceil((imgsiz-1)/2)+1);
img_tmp(1:imgsiz(1),1:imgsiz(2),1:imgsiz(3))=img;
imgsiz_tmp = size(img_tmp);

% get the position of the maximum value of img
if nargin<3 | isempty(MaxValPos)
    [MaxVal,MaxValPosInd] = max(img_tmp(:));     
    [MaxValPos(1),MaxValPos(2),MaxValPos(3)] = ind2sub(imgsiz_tmp,MaxValPosInd); 
end

% shift to center
numDims = ndims(img_tmp);
idx = cell(1, numDims);
for k = 1:numDims
    m = size(img_tmp, k);
    p = ceil((m+1)/2);
    i0 = 1:m;
    i1 = i0 + (MaxValPos(k)-p);
    idi = (i1>=1) & (i1<=m);
    id0{k} = i0(idi);
    id1{k} = i1(idi);
end
img_ctr = zeros(imgsiz_tmp);
img_ctr(id0{:}) = img_tmp(id1{:});

% creat mirror array ...
img_sym = img_ctr;                             
for k = 1:numDims
    img_mir = flipdim(img_sym,k);
    img_sym = max(img_sym,img_mir); 
end

% make symmetric
imgout = img_sym(1:imgsiz(1),1:imgsiz(2),1:imgsiz(3));
