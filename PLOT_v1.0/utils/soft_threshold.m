function x = soft_threshold(y, tau, pentyp, delta, siz)
%
% soft-threshoding for various penalties
%
% gbwang@ucdavis.edu, 12-01-2012
%

switch pentyp
    
    case {'l2','quad'}
        x = y./(1+tau);
        
    case {'l1','hyper'}
        if delta>0
            error('delta cannot be greater than 0');
        end
        x = sign(y) .* max(0, abs(y)-tau);
        
    case {'fair','lange'}
        z = abs(y) - delta - tau;
        x = sign(y) .* ( z + sqrt(z.^2+4*delta*abs(y)) )/2;
        
    case {'tv'}

        if nargin<5 | isempty(siz)
            siz = size(y(:)');
        end            
        y = reshape(y,siz);
        z = sqrt(sum(y.^2,2));
        z = repmat(z,[1 siz(2)]);
        x = y./z.*max(0,z-reshape(tau,siz));
        x(z==0) = 0;
        x = x(:);
end