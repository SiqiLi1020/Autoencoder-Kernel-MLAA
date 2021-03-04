function [x, out, Phi] = eml_dem(yi, ni, G, Gopt, x0, ri, maxit, beta, R)
%--------------------------------------------------------------------------
% Optimization transfer (OT) algorithm (based on expectation maximization)  
% for penalized likelihood image reconstruction of emission data.
% Regularization can be either pixel-based or patch-based. The details of 
% the patch-based regularization and the OT algorithm are described in
%
%   G. Wang and J. Qi, Penalized Likelihood PET Image Reconstruction Using 
%   Patch-Based Edge-Preserving Regularization, IEEE Transactions on Medical
%   Imaging, 2012, 31(12): 2194-2204
%
% This OT algorithm is very similar to De Pierro's MAP-EM (IEEE-TMI 1995) 
% but is extended for patch-based regularization
%--------------------------------------------------------------------------
% INPUT:
%   yi      sinogram in vector
%   ni      multiplicative factor (normalization, attenuation), can be empty
%   G       system matrix
%   Gopt    option set for G, can be empty if G is a matlab sparse matrix
%   ri      addtive factor (randoms and scatters), can be empty
%   x0      initial image estimate, can be empty
%   maxit   maximum iteration number
%   beta    regularization parameter
%   R       regularization operator
%
% OUTPUT
%   x       image estimate in vector
%   out     output
%   Phi     objective function value
%
%--------------------------------------------------------------------------
% Programmer: Guobao Wang @ UC Davis, Qi Lab, 02-01-2012
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
if isempty(maxit)
    maxit = 10;
end

% regularization
if isempty(beta) | nargin<8
    beta = 0;
end
if beta>0
    C1 = sum(abs(R.C),2);
end

% initialization
x    = max(mean(x0(:))*1e-9,x0(:)); x(~Gopt.mask) = 0;
yeps = mean(yi(:))*1e-9;
wx   = Gopt.sens;

% output
if nargin>1
    out = []; Phi = [];
end
out.xest = zeros(length(x(:)), min(maxit,ceil(maxit/Gopt.savestep)+1));
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
    
    % EM update
    yb = ni.*proj_forw(G, Gopt, x) + ri;
    yy = yi./(yb+yeps);
    yy(yb==0&yi==0) = 1;
    xb = proj_back(G, Gopt, ni.*yy);
    xx_em = x ./ wx .* xb;
    xx_em(~Gopt.mask) = 0;
    
    % regularization
    if beta>0
        cx     = R.C*x(:);
        [U, W] = halfQuad(R, R, x, cx);
        gx_reg = R.C'*(W.*cx);           
        wx_reg = abs(R.C)'*(W.*C1);
        xx_reg  = x - gx_reg./wx_reg;
        xx_reg(wx_reg==0) = 0;
    else
        W = []; U = 0; wx_reg = 0; xx_reg = 0;
    end
    
    % objective function value
    if nargout>2
        iy = yb>0;
        Phi(it) = sum(yi(iy).*log(yb(iy))-yb(iy)) - beta*U;
        if it>1 & Phi(it)<Phi(it-1)
            warning('Objective function is not increasing')
        end
    end
    
    % fusion for OT update
    x = eml_prox_sepl2(wx, xx_em, beta*wx_reg, xx_reg);
     
end
if nargout>2 & Gopt.savestep>1
    Phi = Phi([1 Gopt.savestep:Gopt.savestep:end]);
end
proj_clear(Gopt);




    
