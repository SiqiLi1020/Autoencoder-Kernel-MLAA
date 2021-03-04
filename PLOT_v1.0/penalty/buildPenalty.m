function R = buildPenalty(imgsiz, pentyp, penpar, patchsiz, isotropic, R)
%
% Build the regularizer for penalized likelihood image reconstruction
%
% Guobao Wang @ UC Davis (10-01-2012)
%

%% check inputs
imgdim = length(find(imgsiz>0));
if imgdim==1
    imgsiz = [imgsiz(imgsiz(:)>1) 1 1];
elseif imgdim==2
    imgsiz = [imgsiz(imgsiz(:)>0) 1];
end
if nargin<2 | isempty(pentyp)
    pentyp = 'quad';
end
if nargin<3 | isempty(penpar)
    penpar = 1;
end
if nargin<4 | isempty(patchsiz)
    pathsiz = 1;
end
if nargin<5 | isempty(isotropic)
    isotropic = 1;
end
if strcmp(pentyp,'tv')
    if isotropic~=1
        error('tv must set isotropic=1');
    end
end
if nargin<6 | isempty(R)
    R.imgsiz = imgsiz;
end

%% penalty
[potfun, dpot_1, dpot_h] = penalty(pentyp);

%% patch
wlen = 2*floor(patchsiz/2); 
widx = -wlen/2:wlen/2; 
xidx = widx; 
if imgsiz(2)>1
    yidx = widx;
else
    yidx = 0;
end
if imgsiz(3)>1
    zidx = widx;
else
    zidx = 0;
end
n = 1; 
for i = 1:length(xidx)
    for j = 1:length(yidx)
        for k = 1:length(zidx)
            d(n) = xidx(i)^2 + yidx(j)^2 + zidx(k)^2;
            s(:,n) = [xidx(i) yidx(j) zidx(k)]';
            n = n + 1;
        end
    end
end
s = s(:,1:(n-1));
if length(d)==1 d = 1; end
w = 1./sqrt(d); w(d==0) = 2;

J = [1:prod(imgsiz)]';
P = zeros(length(J),size(s,2));
for l = 1:size(s,2)
    P(:,l) = setBoundary3(J, s(:,l), imgsiz);
end

%% output
R.penfun  = potfun;
R.deriv1  = dpot_1;
R.derivh  = dpot_h;
R.pentyp  = pentyp;
R.penpar  = penpar;
R.pshift  = s;
R.pweight = w(:)/sum(w);
R.pindex  = P;
R.isotropic = isotropic;
if ~isfield(R,'spatwgt')
    R.spatwgt = ones(prod(imgsiz),1);
end
R.imgsiz = imgsiz;

%--------------------------------------------------------------------------
function [potfun, dpot_1, dpot_h] = penalty(pentyp)
%--------------------------------------------------------------------------
%
% INPUT 
%   pentyp      type of penalty function: 'quad','huber','lange',...
% OUTPUT
%   potfun      potential function p(t;delta)
%   dpot_1      first derivative of p(t;delta)
%   dpot_2      second derivatie of p(t;delta)    
%   dpot_h      half-quadratic derivative of p(t;delta)

switch lower(pentyp)
    
    case {'quad','gaussian'}
        potfun = inline('w.*(t.^2)/2','w','t','d');
        dpot_1 = inline('w.*t','w','t','d');
        dpot_2 = inline('w','w','t','d');
        dpot_h = inline('w','w','t','d');
        
    case 'huber'
        potfun = @(w,t,d) huber(t, d, 'potfun', w);
        dpot_1 = @(w,t,d) huber(t, d, 'deriv1', w);
        dpot_h = @(w,t,d) huber(t, d, 'derivh', w);
        
    case {'lange','fair'}
        potfun = @(w,t,d) fair(t, d, 'potfun', w);
        dpot_1 = @(w,t,d) fair(t, d, 'deriv1', w);
        dpot_h = @(w,t,d) fair(t, d, 'derivh', w);
    
    case {'hyper','tv'}
        potfun = inline('w.*sqrt(t.^2+d.^2)','w','t','d');
        dpot_1 = inline('w.*t./sqrt(t.^2+d.^2)','w','t','d');
        dpot_2 = inline('w.*d^2./sqrt(t.^2+d^2).^3','w','t','d');
        dpot_h = inline('w./sqrt(t.^2+d^2)','w','t','d');
        
end

%--------------------------------------------------------------------------
function J = setBoundary3(J, d, imgsiz)
%--------------------------------------------------------------------------
[x,y,z] = ind2sub(imgsiz,J);
x = setBoundary1(x+d(1), imgsiz(1)); 
y = setBoundary1(y+d(2), imgsiz(2)); 
z = setBoundary1(z+d(3), imgsiz(3));
J = sub2ind(imgsiz,x,y,z);

%--------------------------------------------------------------------------
function x = setBoundary1(x, N)
%--------------------------------------------------------------------------
if N==1
    x = x.^0;
else
    idx = x(:)>N;
    x(idx) = N - (x(idx)-N);
    idx = x(:)<1;
    x(idx) = 1 + (1-x(idx));
    x = x(:);
end

