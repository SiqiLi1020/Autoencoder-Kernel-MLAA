function [U, W] = halfQuad(H, R, x, cx, penpar)
%--------------------------------------------------------------------------
% Build the optimization transfer weight for pixel-based or patch-based 
% regularizers. 
%
% Guobao Wang @ UC Davis (10-01-2012)
% Last Modified: 03-02-2013
%

% check
if nargin<4
    cx = [];
end
if nargin<5 | isempty(penpar)
    penpar = R.penpar;
end

% if neighborhood is empty
if isempty(H.N)
    U = sum(R.spatwgt.*R.penfun(1,abs(x),penpar));
    W = R.spatwgt.*R.derivh(1,abs(x),penpar); 
    return;
end

% check pixel-based or patch-based
if length(R.pweight)>1
    isPatch = 1;
else
    isPatch = 0;
end

% neighborhood
numpix = prod(R.imgsiz);
J = [1:numpix]';
switch R.isotropic 
    case 1  % Isotropic
    
        % pixel-wise weights
        if isPatch
            d = zeros(numpix,1);
            for k = 1:size(H.N,2)
                d  = d + H.W(:,k).*( ( x(R.pindex(J,:)) - x(R.pindex(H.N(:,k),:)) ).^2 * R.pweight );
            end
        else
            if isempty(cx)
                cx = H.C*x;
            end
            d  = sum(reshape(H.W(:).*cx.^2, [numpix size(H.N,2)]),2);
        end
        d = sqrt(d);
        w = R.spatwgt.*R.derivh(1,d,penpar); 
        U = sum(R.spatwgt.*R.penfun(1,d,penpar));
        
        % weights of neighboring pixels
        if nargout>1
            if isPatch
                W = zeros(numpix,size(H.N,2));
                wk = zeros(size(w));
                for k = 1:size(H.N,2)
                    wk = wk.*0;
                    for l = 1:length(R.pweight)
                        Jl = R.pindex(J,l);
                        wk = wk + w(Jl).*H.W(Jl,k)*R.pweight(l); 
                    end
                    W(:,k) = wk;
                end
            else
                W = H.W.*repmat(w,[1 size(H.N,2)]);
            end
        end
        
    case 0 % Anisotropic
    
        % weights of neighboring pixels
        if isPatch
            if nargout>1
                W = zeros(numpix,size(H.N,2));
            end
            U = 0;
            wk = zeros(numpix,1);
            for k = 1:size(H.N,2)

                % objective function
                d = ( x(R.pindex(J,:)) - x(R.pindex(H.N(:,k),:)) ).^2 * R.pweight;
                U = U + sum(R.spatwgt.*R.penfun(H.W(J,k),sqrt(d),penpar));   

                % weights
                if nargout>1
                    wk = wk.*0;
                    for l = 1:length(R.pweight)
                        Jl = R.pindex(J,l);
                        d  = (x(R.pindex(Jl,:))-x(R.pindex(R.pindex(H.N(:,k),l),:))).^2*R.pweight;
                        wk = wk + R.derivh(H.W(Jl,k),sqrt(d),penpar)*R.pweight(l); 
                    end
                    W(:,k) = R.spatwgt.*wk;
                end
            end
        else
            if isempty(cx)
                cx = H.C*x;
            end
            U  = sum(vec(repmat(R.spatwgt,[1 size(H.W,2)])).*R.penfun(H.W(:), cx, penpar)); 
            if nargout>1
                W = R.derivh(H.W(:), cx, penpar);
            end
        end
end

if nargout>1
    W = W(:);
end
