%--------------------------------------------------------------------------
function put_data(filename, data, type, offset)
%--------------------------------------------------------------------------
% write data into file
%
% gbwang@ucdavis.edu (01-09-2013)
%

if nargin<3 | isempty(type)
    type = 'float32';
end
if nargin<4 | isempty(offset)
    offset = 0;
end

switch type
    case 'int16'
        len = 2;
    case 'float32'
        len = 4;
    case 'int32'
        len = 4;
end

fid = fopen(filename, 'r+');
if fid==-1
    fid = fopen(filename, 'w+');
end
fseek(fid,offset*len,'bof');
fwrite(fid, data, type);
fclose(fid); 