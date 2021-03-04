function LAC = CT2LAC(CT, kVp, method)
% Convert CT HU to linear attenuation coefficients at 511 keV
% 
% Reference: M. Abella et al, Phys Med Biol 2012, 57(9): 2477-2490.
%
% gbwang@ucdavis.edu 10-24-2019
%

if nargin<2 | isempty(kVp)
    kVp = '120';
end
if nargin<3 | isempty(method)
    method = 'trilinear';
end

switch kVp
    case '80'
        a(1) = 9.3e-5;
        b(1) = 0.093;
        a(2) = 3.28e-5;
        b(2) = 0.093;
        a(3) = 4.1e-6;
        b(3) = 0.122;
    case '100'
        a(1) = 9.3e-5;
        b(1) = 0.093;
        a(2) = 4e-5;
        b(2) = 0.093;
        a(3) = 5e-6;
        b(3) = 0.128;
    case '120'
        a(1) = 9.3e-5;
        b(1) = 0.093;
        a(2) = 4.71e-5;
        b(2) = 0.093;
        a(3) = 5.89e-6;
        b(3) = 0.134;
    case '140'
        a(1) = 9.3e-5;
        b(1) = 0.093;
        a(2) = 5.59e-5;
        b(2) = 0.093;
        a(3) = 6.98e-6;
        b(3) = 0.142;
end

% LAC
LAC = zeros(size(CT));

idx = CT>=-1000 & CT<0;
LAC(idx) = a(1)*CT(idx)+b(1);

switch method
    case 'bilinear'
        idx = CT>=0;
        LAC(idx) = a(2)*CT(idx)+b(2);
        
    case 'trilinear'
        idx = CT>=0 & CT<1000;
        LAC(idx) = a(2)*CT(idx)+b(2);
        idx = CT>=1000 & CT<3000;
        LAC(idx) = a(3)*CT(idx)+b(3);
end

