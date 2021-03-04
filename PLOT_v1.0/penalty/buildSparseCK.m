function [C, K] = buildSparseCK(N, W, J, numvox)
%--------------------------------------------------------------------------
% generate sparse matrix C for regularizer. A sparse transform matrix K can
% also be created
%
% Guobao Wang @ UC Davis (10-01-2012)
%

% check input variables
if nargin<3 | isempty(J)
    J = [1:size(N,1)]';
end
if nargin<4 | isempty(numvox)
    numvox = size(N,1);
end

% create C
C = []; 
for n = 1:size(N,2)
    
    % nonzeros
    i = [J; J];
    j = [J; J(N(:,n))];
    s = [ones(numvox,1); -ones(numvox,1)];

    % the sparse matritrix C
    Cn = sparse(i(:), j(:), s(:), numvox, numvox);
    C = [C; Cn];
end

% create the sparse weight matrix K
if nargout>1
    i = repmat(J,[1 size(N,2)]);
    j = J(N);
    K = sparse(J(i(:)), J(j(:)), W(:), numvox, numvox);
end
