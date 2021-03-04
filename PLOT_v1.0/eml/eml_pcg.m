function [x, out, Phi] = eml_pcg(yi, ni, G, Gopt, x0, ri, maxit, beta, R)
%--------------------------------------------------------------------------
% Preconditioned conjugate gradient algorithm for penalized likelihood
% image reconstruction of emission data. References:
%
% 1. L. Kaufman, Implementing and accelerating the EM algorithm for
%    positron emission tomography. IEEE Transations on Medical Imaging,
%    6(1):37-51, 1987.
% 2. E. U. Mumcuoglu, R. Leahy, S. R. Cherry, and Z. Y. Zhou, Fast
%    gradient-based methods for Bayesian reconstruction of transmission
%    and emission PET images, IEEE Transactions on Medical Imaging, 13(4): 
%    687-701, 1994.
% 3. J. A. Fessler, Conjugate-gradient preconditioning methods for 
%    shift-variant PET image reconstruction, IEEE Transactions on Image
%    Processing, 8(5): 688-699, 1999 
%
%--------------------------------------------------------------------------
% INPUT:
%   yi      sinogram in vector
%   ni      multiplicative factor (normalization, attenuation), can be empty
%   G       system matrix
%   Gopt    option set for G, can be empty if G is a matlab sparse matrix
%   x0      initial image estimate, can be empty
%   ri      addtive factor (randoms and scatters), can be empty
%   maxit   maximum iteration number
%   beta    regularization parameter
%   R       regularization operator
%
% OUTPUT
%   x       image estimate in vector
%   out     output
%   Phi     objective function
%
%--------------------------------------------------------------------------
% Programmer: Guobao Wang @ UC Davis, 12-01-2012
% Last Modified: 06-28-2013
%--------------------------------------------------------------------------

%% check inputs
imgsiz = R.imgsiz;
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end
[yi, ri, ni] = sino_preprocess(yi, ri, ni);

% set Gopt
Gopt = setGopt(ni, G, Gopt);
if ~isfield(Gopt,'precond') | isempty(Gopt.precond)
    Gopt.precond = 'em';
end
if ~isfield(Gopt,'conj') | isempty(Gopt.conj)
    Gopt.conj = 1;
end
if isempty(maxit)
    maxit = 10;
end

% regularization
if isempty(beta) | nargin<8
    beta = 0;
end
R = setRopt(R);
if beta>0
    C1 = sum(abs(R.C),2);
end

% initialization
x    = max(mean(x0(:))*1e-9,x0(:)); x(~Gopt.mask) = 0;
yb   = ni.*proj_forw(G, Gopt, x) + ri;
cx   = R.C*x;
yeps = mean(yi(:))*1e-9;
wx   = Gopt.sens;

% output
nnnum = 0;
rdnum = 0;
step  = 1;
if nargin>1
    out = []; Phi = [];
end
out.xest = zeros(length(x(:)),  min(maxit,ceil(maxit/Gopt.savestep)+1));
t1 = tic;

%% iterative loop
for it = 1:maxit     
    
    % save data
    if Gopt.disp
        disp(sprintf('iteration %d',it));
    end
    if nargout>1 & ( it==1 | rem(it,Gopt.savestep)==0 )
        itt = min(it, floor(it/Gopt.savestep) + 1);
        out.step(itt)   = it;
        out.time(:,itt) = toc(t1);        
        out.xest(:,itt) = x;
    end
    
    % gradient
    yy = yi./(yb+yeps);
    yy(yb==0&yi==0) = 1;
    gx = proj_back(G, Gopt, ni.*(yy-1));
    if beta>0
        [U, W] = halfQuad(R, R, x, cx);
        gx_reg = R.C'*(W.*cx);           
        gx     = gx - beta*gx_reg;
    else
        W = []; U = 0; cx = 0; dc = 0;
    end
    gx(~Gopt.mask) = 0;

    % objective function value
    if nargout>2
        iy = yb>0;
        Phi(it) = sum(yi(iy).*log(yb(iy))-yb(iy));
        if beta>0
            Phi(it) = Phi(it) - beta*U;
        end
        if it>1 & Phi(it)<Phi(it-1)
            warning(sprintf('Objective function is not increasing, %2.1g, step=%2.1e',Phi(it)-Phi(it-1),step))
        end    
    end
    
    % preconditioned direction    
    switch Gopt.precond
        case 'no'
            px = gx;
                       
        case 'em'
            px = x./wx.*gx;
            
        otherwise
            error('unknown preconditioner');
    end
    px(~Gopt.mask) = 0;
    
    % conjugate direction
    if ( it==1 | Gopt.conj==0 ) 
        dx = px;
    else
        gamma_den = sum(gx_old.*px_old,1);
        gamma = sum((gx-gx_old).*px,1) ./ gamma_den;
        if gamma_den==0 gamma=0; end;
        dx = px + gamma*dx;
    end
    gx_old = gx;
    px_old = px;
    
    % truncated direction
    if Gopt.nonneg
        dx = truncate_direction(x, dx);
    end
    
    % reset direction if needed
    cosdxgx = (dx'*gx)/norm(dx,2)/norm(gx,2);
    if cosdxgx<0.001
        if Gopt.debug
            disp(sprintf('cosdxgx=%2.1e, resetting direction',cosdxgx));
        end
        rdnum = rdnum + 1;
        dx = gx;
    end
    
    % increments
    dy = ni.*proj_forw(G, Gopt, dx);
    if beta>0
        dc = R.C * dx;
    end
    
    % step size by line search
    step = eml_line_search(yi, yb, dy, x, dx, beta, W, cx, dc, inf, R);
    
    % nonnegativity
    if Gopt.nonneg
        x_new = x + step*dx;
        idx = x_new < 0;
        if any(idx)
            if Gopt.debug
                temp = 100*length(find(idx))/length(find(Gopt.mask));
                disp(sprintf('image has %3.1f%% negatives, running a 2nd line search',temp));
            end
            nnnum = nnnum + 1;
            x_new(idx) = 0;
            dx = x_new - x;
            dy = ni.*proj_forw(G, Gopt, dx);
            if beta>0
                dc = R.C * dx;
            end
            step = eml_line_search(yi, yb, dy, x, dx, beta, W, cx, dc, 1.0, R);
        end
    end
    
    % update
    x  = x  + step*dx;
    yb = yb + step*dy;
    cx = cx + step*dc;
    
end

% display
if Gopt.debug
    disp(sprintf('Total number of second line search  : %d', nnnum));
    disp(sprintf('Total number of resetting direction : %d', rdnum));
end
if nargout>2 & Gopt.savestep>1
    Phi = Phi([1 Gopt.savestep:Gopt.savestep:end]);
end
proj_clear(Gopt);

    
