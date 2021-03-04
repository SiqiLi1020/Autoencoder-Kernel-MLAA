function R = setRopt(R)
% set options for R
%

% for line search
if ~isfield(R,'lstype') | isempty(R.lstype)
    R.lstype = 'orig';
end
if ~isfield(R,'lsiter') | isempty(R.lsiter)
    R.lsiter = 5;
end

