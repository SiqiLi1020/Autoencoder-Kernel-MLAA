function [x, out, Phi] = eml_tot(yi, ni, G, Gopt, x0, ri, maxit, beta, R, sigma)
%--------------------------------------------------------------------------
% A trust optimization transfer algorithm for penalized likelihood PET 
% image reconstruction using edge-preserving regularization. 
%
%   G. Wang and J. Qi, Penalized likelihood edge-preserving PET image
%   reconstruction using trust optimization transfer, 12th International 
%   Meeting on Fully Three-Dimensional Image Reconstruction in Radiology 
%   and Nuclear Medicine (Fully3D), Lake Tahoe, June 16-21 2013, pp. 70-73.
%
%--------------------------------------------------------------------------
% INPUT:
%   yi      sinogram
%   ni      multiplicative factor (normalization, attenuation), can be empty
%   G       system matrix
%   Gopt    option set for G, can be empty if G is a matlab sparse matrix
%   x0      initial image estimate, can be empty
%   ri      addtive factor (randoms and scatters), can be empty
%   maxit   maximum iteration number
%   beta    regularization parameter
%   R       regularization operator
%   sigma   inital sigma value for the trust surrogate
%
% OUTPUT
%   x       image estimate
%   out     output
%   Phi     objective function
%
%--------------------------------------------------------------------------
% Programmer: Guobao Wang @ UC Davis, 12-11-2012
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
    Gopt.precond = 'ot';
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
if ~isfield(R,'lstype') | isempty(R.lstype)
    R.lstype = 'surr';
end
R = setRopt(R);
C1 = sum(abs(R.C),2);
if nargin<10
    sigma = [];
end

% get initial sigma
if isempty(sigma)
    x = proj_back(G, Gopt, ni.*(yi-ri))./Gopt.sens;
    x(x<mean(x(x>0))) = 0;
    x(~Gopt.mask) = 0;
    yy = ni.*proj_forw(G, Gopt, x);
    sc = sum(yy.*(yi-ri))/sum(yy.^2);
    sigma = 0.1 * (sc*mean(x(x>0)));    
end
disp(sprintf('Initial sigma in TrustOT is %3.2g. \n User can tune this parameter for faster convergence.',sigma));

% initialization
x    = max(mean(x0(:))*1e-9,x0(:)); x(~Gopt.mask) = 0;
yb   = ni.*proj_forw(G, Gopt, x) + ri;
cx   = R.C*x;
yeps = mean(yi(:))*1e-9;
wx   = Gopt.sens;

% setting for monotonic continuation
deltaOrig = R.penpar;
sigma     = max(sigma, deltaOrig);
R.penpar  = sigma;
rho       = 1;
restartIt = 1;

% initial objective function
iy = yb>0; 
L  = sum(yi(iy).*log(yb(iy)))-sum(yb(iy));
U  = halfQuad(R, R, x, cx, deltaOrig);

% output
nnnum = 0;
rdnum = 0;
if nargin>1
    out = []; Phi = [];
end
out.xest = zeros(length(x(:)), min(maxit,ceil(maxit/Gopt.savestep)+1));
t1 = tic; 

%% iterative loop
for it = 1:maxit     
    
    % save data
    if Gopt.disp & ~Gopt.debug
        disp(sprintf('iteration %d',it));
    end
    if nargout>1 & ( it==1 | rem(it,Gopt.savestep)==0 )
        itt = min(it, floor(it/Gopt.savestep) + 1);
        out.step(itt) = it;
        out.time(:,itt) = toc(t1);
        out.xest(:,itt) = x;
        out.sigma(itt)  = sigma;
    end
    
    % EM update
    if rho>0
        yy = yi./(yb+yeps);
        yy(yb==0&yi==0) = 1;
        xb = proj_back(G, Gopt, ni.*yy);
        xx_em = x ./ wx .* xb;
        xx_em(~Gopt.mask) = 0;
    end
    
    % regularization
    [Q, W] = halfQuad(R, R, x, cx);
    gx_reg = R.C'*(W.*cx);
    wx_reg = abs(R.C)'*(W.*C1);
    xx_reg = x - gx_reg./wx_reg;
        
    % fusion for OT update
    if strcmp(Gopt.precond,'ot')
        xx = eml_prox_sepl2(wx, xx_em, beta*wx_reg, xx_reg);
    end
    
    % objective function value
    Phi(it) = L - beta*U;
    F       = L - beta*Q;
    if Gopt.disp & it>1 & Phi(it)<Phi(it-1)
        warning(sprintf('Objective function is not increasing, %2.1g',Phi(it)-Phi(it-1)));    
    end
    
    % gradient
    gx = ( xb - wx ) - beta*gx_reg;
    gx(~Gopt.mask) = 0;
    
    % (implicitly or explicitly) preconditioned direction    
    switch Gopt.precond
        case 'no'
            px = gx;
            
        case 'ot'
            px = xx - x; 
                       
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
    cosdxgx = (dx'*px)/norm(dx,2)/norm(px,2);
    if cosdxgx<0.001
        if Gopt.debug
            disp(sprintf('cosdxgx=%2.1e, resetting direction',cosdxgx));
        end
        rdnum = rdnum + 1;
        if strcmp(Gopt.precond,'ot')
            dx = px;
        else
            dx = gx;
        end
        dx = gx;   % useful when delta islarge (e.g. delta=1)
    end
    
    % increments
    dy = ni.*proj_forw(G, Gopt, dx);
    dc = R.C * dx;
    
    % step size by line search
    step = eml_line_search(yi, yb, dy, x, dx, beta, W, cx, dc, inf, R);
    if strcmp(Gopt.precond,'ot') & (Gopt.conj==0)
        step = max(1,step);
    end
    
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
            dc = R.C * dx;
            step = eml_line_search(yi, yb, dy, x, dx, beta, W, cx, dc, 1.0, R);
        end
    end
    
    % gain factor
    x_new = x  + step*dx;
    y_new = yb + step*dy;
    c_new = cx + step*dc;
    iy    = y_new>0;
    L_new = sum(yi(iy).*log(y_new(iy)))-sum(y_new(iy));
    Q_new = halfQuad(R, R, x_new, c_new);
    F_new = L_new - beta*Q_new;
    U_new = halfQuad(R, R, x_new, c_new, deltaOrig);
    Phi_new = L_new - beta*U_new;
    rho   = (Phi_new-Phi(it))/(F_new-F);   
    nuv   = rho*(Phi_new-Phi(it))/(Phi_new-Phi(restartIt));
    if Gopt.debug
        disp(sprintf('iter %d, sigma: %2.1e, rho=%2.1f, restart at: %d, nuv=%2.1e',it, sigma, rho,restartIt,nuv));
    end
    
    % update
    if rho>0 & (F_new>F)
        x  = x_new;  
        yb = y_new; 
        cx = c_new;
        L  = L_new;
        U  = U_new;
    end
    
    % update the dampling parameter
    if rho<0 | (F_new<F)
        sigma = max(deltaOrig, sigma/3);
        restartIt = it;
        if Gopt.debug & sigma>deltaOrig
            disp(sprintf('step=%2.1e, dF=%2.1e, test a new sigma = %2.1e',step, F_new-F, sigma)) 
        end
    elseif nuv<1e-2
        sigma = max(deltaOrig, sigma/3);
        restartIt = it+1;
    end
    R.penpar = sigma;
        
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




    
