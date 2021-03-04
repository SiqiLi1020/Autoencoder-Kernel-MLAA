%--------------------------------------------------------------------------
function [data, datalen] = get_data(filename, size, type, offset)
%--------------------------------------------------------------------------
% read data from file
%
% gbwang@ucdavis.edu (01-09-2013)
%

if nargin<2 | isempty(size)
    size = inf;
end
if nargin<3 | isempty(type)
    type = 'float32';
end
if nargin<4 | isempty(offset)
    offset = 0;
end
if isempty(filename)
    data = [];
    datalen = 0;
    return
end

% length of element
switch type
    case 'int16'
        len = 2;
    case 'float32'
        len = 4;
    case 'int32'
        len = 4;
end

fid = fopen(filename, 'r');
if fid==-1
    disp(sprintf('WARNING: There Is No File %s', filename));
    data = [];
    datalen = 0;
else
    % pointer
    fseek(fid, offset*len, 'bof');
        
    [data, datalen] = fread(fid, size, type);    
    fclose(fid);
end
