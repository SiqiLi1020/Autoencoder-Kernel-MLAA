 function f = huber(t, d, type, w)
%--------------------------------------------------------------------------
% Huber potential function
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if nargin<3 | isempty(type)
    type = 'potfun';
end
switch type
    case 'potfun'
        f = t.^2 / 2;
        i = abs(t) > d;
        f(i) = d * abs(t(i)) - d.^2/2;
        
    case 'deriv1'
        f = t;
        i = abs(t) > d;
        f(ii) = d * sign(t(i));
        
    case 'derivh'
        f = ones(size(t));
        i = abs(t) > d;
        f(i) = d ./ abs(t(i));

    otherwise
        error('unknown type for Hunber potential')
end

if nargin>3
    f = f.*w;
end
