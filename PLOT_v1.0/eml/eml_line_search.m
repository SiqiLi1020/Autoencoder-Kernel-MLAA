function step = eml_line_search(yi, yb, dy, x, dx, alpha, w, cx, dc, maxstep, R)
%--------------------------------------------------------------------------
% Line search to find the step size in penalized likelihood reconstruction
% of emission data.
%
% Guobao Wang @ UCDavis, 12-01-2011
% Last Modified: 03-03-2013
%

% check
if nargin<11 | ~isfield(R,'lstype') | isempty(R.lstype) 
    R.lstype = 'surr';
end
if nargin<11 | ~isfield(R,'lsiter') | isempty(R.lsiter)
    R.lsiter = 5;
end

% must use quadratic surrogate for patch-based regularization
if ~isfield(R,'pweight') | length(R.pweight)>1
    R.lstype = 'surr';
end

% determine the upper limit of the step size
idy = find(dy<0);
if ~isempty(idy)
    tempmax = - max(yb(idy)./dy(idy));    
else
    tempmax = inf;    
end
maxstep = min(tempmax, maxstep);

% initilization
yeps = mean(yi)*1e-9;
if alpha>0
    if isempty(w)
        [U, w] = halfQuad(R, R, x, cx);
    end
    dcwdc = dc'*(w.*dc);
    cxwdc = cx'*(w.*dc);
end
mu   = 1;
step = 0;

% iteration loop
for it = 1:R.lsiter
    
    % derivatives of likelihood function
    yy  = yb + step * dy;
    temp= 1-yi./max(yy,yeps);
    iy  = yi==0 & yy==0; temp(iy) = 0; % This is very important to get monotonicity
    g1  = sum(dy.*temp,1);
    temp= yi./max(yy,yeps).^2; temp(iy) = 1;
    g2  = sum(dy.^2.*temp,1);
    
    % derivatives of spatial penalty funtion
    if alpha>0
        g1 = g1 + alpha * ( cxwdc + step*dcwdc ); 
        g2 = g2 + alpha * dcwdc;
    end
    
    % Newton-Raphson update
    step = step - mu .* g1./g2;
    if g2==0 step = 0; end
    if step>maxstep;
        step = maxstep;
        mu = mu/2;
    elseif step<0
        step = step*0.01;
        mu = mu/2;
    end
    
    % update
    if alpha>0 & strcmp(R.lstype,'orig')
        [U, w] = halfQuad(R, R, x+step*dx, cx+step*dc);
        dcwdc = dc'*(w.*dc);
        cxwdc = cx'*(w.*dc);
    end
    
end
