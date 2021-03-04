function Phi = eml_objfun(yi, ni, G, Gopt, xs, ri, alpha, R)
%
% calculate the cost function value of the penalized maximum likelihood
%
% Guobao Wang @ UC Davis, 12-01-2012
%

%% check inputs
imgsiz = R.imgsiz;
numpix = prod(imgsiz);

if isempty(Gopt) | ~isfield(Gopt,'mtype')
    Gopt.mtype = 'sparse';
end
if isempty(G) & strcmp(Gopt.mtype,'sparse')
    error('System matrix G must be provided')
end
if ~isfield(Gopt,'mask') | isempty(Gopt.mask)
    Gopt.mask = ones(numpix,1)>0;
end
[yi, ri, ni] = sino_preprocess(yi, ri, ni);

%% iterative loop
for i = 1:size(xs,2)
    x  = max(0,xs(:,i)); 
    yb = ni.*proj_forw(G, Gopt, x) + ri;
    iy = yb>0;
    Phi(i) = sum(yi(iy).*log(yb(iy))-yb(iy));
    if alpha>0
        U = halfQuad(R, R, x);
        Phi(i) = Phi(i) - alpha*U;
    end    
end







    
