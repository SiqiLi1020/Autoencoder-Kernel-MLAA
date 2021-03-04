 function f = fair(t, d, type, w)
%--------------------------------------------------------------------------
% Fair (a.k.a. Lange) potential function
%
% Guobao Wang @ UC Davis (02-01-2013)
%

if nargin<3 | isempty(type)
    type = 'potfun';
end
switch type
    case 'potfun'
        if d==0
            f = abs(t);
        else
            f = d.*(abs(t)/d-log(1+abs(t)/d));
        end
        
    case 'deriv1'
        if d==0
            f = sign(t);
        else
            f = t./(abs(t)+d);
        end
        
    case 'deriv2'
        f = -d./(abs(t)+d).^2;
        
    case 'derivh'
        f = 1./(abs(t)+d);

    otherwise
        error('unknown type for Hunber potential')
end

if nargin>3
    f = f.*w;
end
