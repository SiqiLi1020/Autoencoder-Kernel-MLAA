function [u, x, out] = psct_kmlaa(yi, ni, A, Aopt, u, G, Gopt, x, ri, maxit,K_CT);
%
% CT-kernelized MLAA algorithm for synergistic reconstructon of time-of-flight 
% PET activity and attenuation images
%
% gbwang@ucdavis.edu Aug 02, 2019
%
% Last updated: Aug 11 2020
%

%% check inputs for reconstruction

% total frame number
numprj = prod(Aopt.prjsiz);
numbin = size(yi,1)/numprj;

% normalization
if nargin<3 | isempty(ni)
    ni = ones(size(yi));
end

% preprocess sinograms
[yi, ri, ni] = sino_preprocess(yi, ri, ni);

% CT operator
Aopt = setGopt(ones(numprj,1), A, Aopt);
numpix_CT = prod(Aopt.imgsiz);
if nargin<5 | isempty(u)
    u = zeros(numpix_CT,1);
end

% PET operator
numpix = prod(Gopt.imgsiz);
Gopt = setGopt(ni, G, Gopt);
if nargin<8 | isempty(x)
    x = ones(numpix,1);
end
if Gopt.emisFlag
    x = max(mean(x(:))*1e-9,x(:)); 
    x(~Gopt.mask,:) = 0;
end
if nargin<9 | isempty(ri)
    yeps = mean(yi(:))*1e-9;
    ri = ones(size(yi))*yeps;
    yi = yi + ri;
end

% iteration number
if nargin<10 | isempty(maxit)
    maxit = 10;
end

% kernel for gamma-ray CT reconstruction
if nargin<11 | isempty(K_CT)
    K_CT = speye(numpix_CT);
end
Aopt.kernel = K_CT;

% iteration number
if isfield(Aopt,'MaxIt')
    it_attn = Aopt.MaxIt;
else
    it_attn = 1;
end
if isfield(Gopt,'MaxIt')
    it_emis = Gopt.MaxIt;
else
    it_emis = 1;
end

% initialization
if ~Aopt.attnFlag
    li_attn = proj_forw(A, Aopt, u);
    ai = exp(-li_attn); 
    ai = repmat(ai,[1 numbin]);
    mi = ni.*ai(:);
    gg = proj_back(G, Gopt, mi);
else
    aa = proj_forw(A, Aopt, ones(numpix_CT,1)); % aa = A*1
end
if ~Gopt.emisFlag
    li_emis = proj_forw(G, Gopt, x);
end

% output
if nargin>1
    out = []; L = [];
end
out.xest = zeros(length(x(:)), min(maxit,ceil(maxit/Gopt.savestep)+1));
out.uest = zeros(length(u(:)), min(maxit,ceil(maxit/Gopt.savestep)+1));
t1 = tic;

%% iterative loop
for it = 1:maxit     
        
    % save data
    if nargout>1 & ( it==1 | rem(it,Gopt.savestep)==0 )
        itt = min(it, floor(it/Gopt.savestep) + 1);
        if Gopt.disp==1
            disp(sprintf('iteration %d',it));
        end
        out.step(itt)   = it;
        out.time(:,itt) = toc(t1);        
        out.xest(:,itt) = x(:);
        out.uest(:,itt) = u(:);
    end
	
    % forward projection
    if Gopt.emisFlag & it==1
        li_emis = proj_forw(G, Gopt, x);
    end
    if Aopt.attnFlag & it==1
        li_attn = proj_forw(A, Aopt, u);
        ai = exp(-li_attn); 
        ai = repmat(ai,[1 numbin]);
    end
    
    % likelihood function
    yb = ni.*ai(:).*li_emis + ri;
    L(it) = sum(yi(:).*log(yb(:))-yb(:));
    
    % attenuation estimate
    if Aopt.attnFlag 
        mi = ni.*li_emis;
        [u, V] = attn_tofml_sps(yi, mi, A, Aopt, u, ri, it_attn, li_attn, aa);
        li_attn = proj_forw(A, Aopt, u);
        ai = repmat(exp(-li_attn), [1 numbin]);
    end
    
    % emission estimate
    if Gopt.emisFlag
        if Aopt.attnFlag 
            mi = ni.*ai(:);
            gg = proj_back(G, Gopt, mi);
        end
        [x, U] = emis_ml_em(yi, mi, G, Gopt, x, ri, it_emis, li_emis, gg);   
        li_emis = proj_forw(G, Gopt, x);
    end
        
end

% output
out.cost = L;

proj_clear(Gopt);
