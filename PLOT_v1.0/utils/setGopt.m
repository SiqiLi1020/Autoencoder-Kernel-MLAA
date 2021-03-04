function Gopt = setGopt(ni, G, Gopt)
% set the option for Gopt
%
%

% system matrix
if isempty(Gopt) | ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end 
if isfield(Gopt,'mtype') & strcmp(Gopt.mtype, 'matlab') & isempty(G)
    error('G must be provided for the matlab sparse matrix type')
end

% algorithm
if ~isfield(Gopt,'nonneg') | isempty(Gopt.nonneg)
    Gopt.nonneg = 1;
end

% debug, display, save
if ~isfield(Gopt,'disp') | isempty(Gopt.disp)
    Gopt.disp = 1;
end
if ~isfield(Gopt,'debug') | isempty(Gopt.debug)
    Gopt.debug = 0;
end
if ~isfield(Gopt,'savestep') | isempty(Gopt.savestep)
    Gopt.savestep = 10;
end

% sensitivity map
if ~isfield(Gopt,'sens') | isempty(Gopt.sens)
    Gopt.sens = proj_back(G, Gopt, ni);
end

% mask
if ~isfield(Gopt,'mask') | isempty(Gopt.mask)
    Gopt.mask = Gopt.sens>0;
end
