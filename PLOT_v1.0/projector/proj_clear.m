function proj_clear(Gopt)
%--------------------------------------------------------------------------
% delete all temporal files 
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'sparse';
end

switch Gopt.mtype
 
    case 'usc'
        delete temp.image.usc;
        delete temp.image.usc.forw;
        delete temp.proj.usc;
        delete temp.proj.usc.back;
        
    case 'zhou'
        delete temp.image;
        delete temp.image.forw;
        delete temp.proj;
        delete temp.proj.back;       
end
