function xi = eml_prox_sepl2(wi, yi, bet, zi, ri)
%--------------------------------------------------------------------------
% estimate an intermediate estimate from Poisson data
% xi = argmax wi * ( yi*log(xi+ri) - (xi+ri) ) - bet/2 * (xi-zi)^2;
%
% Guobao Wang @ UC Davis
%

% check
if isempty(wi)
    wi = ones(size(yi));
end
if isscalar(bet)
    bet = bet*ones(size(yi));
end
if nargin<5 | isempty(ri)
    ri = zeros(size(yi));
end

% solve the quadratic equation
bi = wi - bet.*(zi+ri);
di = sqrt(bi.^2+4*bet.*wi.*yi);

xi = zeros(size(yi));
ib = bet==0;
xi(ib) = yi(ib);

ib = bet>0 & bi>0;
xi(ib) = 2*wi(ib).*yi(ib)./(di(ib)+bi(ib));

ib = bet>0 & bi<=0;
xi(ib) = (di(ib)-bi(ib))./(2*bet(ib));

xi = xi - ri;
