function [x, yb, cx, Phi] = pwls_iter(wy, yi, ni, G, Gopt, x0, ri, maxit, alpha, u, C, beta, v, yb, cx)
%
% Penalized weighted least squares (PWLS) image upate for minimizing the 
% cost function defined by
%
%  x = argmin_x 1/2(y-Ax)'W(y-Ax) + alpha/2*||u-Cx||^2 + beta/2*||v-x||^2
%
% where A = Diag[ni]*G and W = Diag[wy]; Preconditioned conjugate gradient 
% is used. Note that this is usually a subroutine in alternating direction. 
%
% gbwang@ucdavis.edu, 12-01-2012
%

%% check inputs
imgsiz = Gopt.imgsiz;
numpix = prod(imgsiz);

if isempty(wy)
    wy = ones(size(yi));
end
if isempty(ni)
    ni = ones(size(yi));
end
if isempty(ri)
    ri = zeros(size(yi));
end
if isempty(Gopt) | ~isfield(Gopt,'mtype')
    Gopt.mtype = 'sparse';
end    
if isempty(G) & strcmp(Gopt.mtype,'sparse')
    error('System matrix G must be provided')
end
if ~isfield(Gopt,'mask') | isempty(Gopt.mask)
    Gopt.mask = ones(numpix,1)>0;
end
if isempty(x0)
    x0 = ones(numpix,1)*0.01;
end
if isempty(maxit)
    maxit = 2;
end
if isempty(alpha)
    alpha = 0;
end
if isempty(u)
    u = zeros(size(x0));
end
if alpha>0 & (nargin<11 |isempty(C))
    error('C must be provided');
end
if ~isfield(Gopt,'precond') | isempty(Gopt.precond)
    Gopt.precond = 'cone';
end
if ~isfield(Gopt,'nonneg') | isempty(Gopt.nonneg)
    Gopt.nonneg = 0;
end
if nargin<12 | isempty(beta)
    beta = 0;
end
if nargin<13 | isempty(v)
    v = zeros(size(x0));
end
if nargin<14
    yb = [];
end
if nargin<15
    cx = [];
end

% preconditioner
if  strcmp(Gopt.precond,'ot') | strcmp(Gopt.precond,'em') 
    if ( ~isfield(Gopt,'wc') | isempty(Gopt.wc) )
        yy = ni.*proj_forw(G, Gopt, ones(numpix,1));
        wx = proj_back(G, Gopt, ni.*wy.*yy);
        c  = sum(abs(C),2);    
        wc = abs(C)'*c;
    else
        wx = Gopt.wx;
        wc = Gopt.wc;
    end
end
if strcmp(Gopt.precond,'cone') 
    if ( ~isfield(Gopt,'lam') | isempty(Gopt.lam) )
        lam = cone_coef(wy, G, Gopt, ni, alpha, [], C, beta);
    else
        lam = Gopt.lam;
    end
    lam = lam + max(lam(:))*0.001;
end
mask = Gopt.mask;

% initialization
x = x0(:); x(~mask) = 0;
if isempty(yb)
    yb = ni.*proj_forw(G, Gopt, x) + ri;
end
if alpha>0 & isempty(cx)
    cx = C*x;
end

%% iterative loop
for it = 1:maxit     
    
    % objective function value
    if nargout>3
        Phi(it) = sum(wy.*(yi-yb).^2)/2;
        if alpha>0
            Phi(it) = Phi(it) + alpha/2*sum((cx-u).^2);
        end
        if beta>0
            Phi(it) = Phi(it) + beta/2*sum((x-v).^2);
        end
        if it>1 & Phi(it)>Phi(it-1)
            error(sprintf('objective function is not decreasing at iter %d',it));
        end
    end
    
    % gradient
    gx = proj_back(G, Gopt, ni.*wy.*(yb-yi));
    if alpha>0
        gx = gx + alpha*( C'*(cx-u) );           
    end    
    if beta>0
        gx = gx + beta*(x-v);
    end
    gx(~mask) = 0;       
    
    % preconditioning
    switch Gopt.precond
        case 'no'
            qx = 1;
            px = -gx;
            
        case 'em'
            qx = wx;
            px = (-gx) ./ qx;
            px(qx==0) = 0;
            
        case 'ot'
            qx = wx;
            if alpha>0
                qx = qx + alpha*wc;
            end
            if beta>0
                qx = qx + beta;
            end
            px = (-gx) ./ qx;
            px(qx==0) = 0;
                       
        case 'cone'
            gx_fft = fftn(reshape(-gx,imgsiz));
            px = ifftn(gx_fft./lam);
            px = real(px(:));
            
        otherwise
            error('unknown preconditioner');
    end
    px(~mask) = 0;
    
    % conjugate direction
    if it==1
        dx = px;
    else
        gamma = sum((gx-gx_old).*px,1) ./ sum(gx_old.*px_old,1);
        dx = px + dx * diag(gamma);
    end
    gx_old = gx;
    px_old = px;
        
    % truncated direction
    if Gopt.nonneg
        dx = truncate_direction(x, dx);
    end

    % reset direction if needed
    if abs(dx'*gx)/norm(dx,2)/norm(gx,2)<0.001
        disp('resetting direction')
        dx = -gx;
    end

    % increments
    dy = ni.*proj_forw(G, Gopt, dx);
    if alpha>0
        dc = C*dx;
    end
    
    % step size
    step_num = sum(-dx.*gx);
    step_den = sum(wy.*dy.^2);
    if alpha>0
        step_den = step_den + alpha*sum(dc.^2);
    end
    if beta>0
        step_den = step_den + beta*sum(dx.^2);
    end
    step = step_num/step_den;
    
    % nonnegativity
    if Gopt.nonneg
        xnew = x + step*dx;
        idx  = xnew < 0;
        if any(idx)
            xnew(idx) = 0;
            dx = xnew - x;
            dy = ni.*proj_forw(G, Gopt, dx);
            step_num = sum(-dx.*gx);
            step_den = sum(wy.*dy.^2);
            if alpha>0
                step_den = step_den + alpha*sum(dc.^2);
            end
            if beta>0
                step_den = step_den + beta*sum(dx.^2);
            end
            step = step_num/step_den;
            step = min(1, step);
        end
    end
    
    % update
    x  = x  + step*dx;  
    yb = yb + step*dy;
    if alpha>0
        cx = cx + step*dc;
    end
    
end








    