function dx = truncate_direction(xt, dx)
%--------------------------------------------------------------------------
% truncate the update direction to ensure nonnegativity on images. x must
% be a column vector or a matrix.
%
% gbwang@ucdavis.edu (01-09-2013)
%

for m = 1:size(xt, 2)
    xm = xt(:,m);
    dm = dx(:,m);
    idx = xm < mean(xm)/1e9;
    idd = dm < 0;
    idx = idx & idd;
    dm(idx) = 0;
    dx(:,m) = dm;
end