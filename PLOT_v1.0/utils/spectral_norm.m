function [s, lam] = spectral_norm(wi, ni, G, Gopt, alpha, wk, C, maxit, beta)
%
% Power method to compute spetral norm ||G'D[wi]G + alpha*C'D[wk]C' + beta||_2
%
% gbwang@ucdavis.edu (01-09-2013)
%

%% check
if nargin<4 | isempty(Gopt) | ~isfield(Gopt,'prjsiz')
    M = size(G,1);
    N = size(G,2); 
else
    M = prod(Gopt.prjsiz);
    N = prod(Gopt.imgsiz);
end

if isempty(wi)
    wi = ones(M,1);
end
if isempty(ni)
    ni = ones(M,1);
end
if nargin<4
    Gopt = [];
end
if nargin<5 | isempty(alpha)
    alpha = 0;
end
if nargin<6 | isempty(wk)
    wk = 1;
end
wk = wk(:);
if nargin<8 | isempty(maxit)
    maxit = 100;
end
if nargin<9 | isempty(beta)
    beta = 0;
end
    
% initial
x = rand(N,1);
f = proj_back(G, Gopt, wi.*ni.^2.*proj_forw(G, Gopt, x));
if alpha>0
    f = f + alpha*(C'*(wk.*(C*x)));
end
if beta>0
    f = f + beta*x;
end

%% iteration
for it = 1:maxit
    
    % normalization
    x = f/norm(f);
    
    % matrix-vector product f=K'Kx
    f = proj_back(G,Gopt,wi.*ni.^2.*proj_forw(G,Gopt,x));
    if alpha>0
        f = f + alpha*(C'*(wk.*(C*x)));
    end
    if beta>0
        f = f + beta*x;
    end
    
    % norm
    lam(it) = sqrt(x'*f);
    
    if ( it>1 & (lam(it)-lam(it-1))/(lam(it)-lam(1))<1e-4 ) | it==maxit
        s = lam(end);
        break;
    end
end

    
    