function [v, d] = image_roi(x, roi, typ, x0)
% extract information from a region of interest
%
% gbwang@ucdavis.edu Jan. 2012
%
if isvector(x)
    x = x(:);
end
if nargin<3 | isempty(typ)
    typ = 'mean';
end

roi = roi>0;
for k = 1:size(roi,2)
    switch typ
        case 'mean'
            v(k) = mean(x(roi(:,k),:));
            d(k) = std(x(roi(:,k),:),0,1);

        case 'max'
            v(k) = max(x(roi(:,k),:));
            d(k) = std(x(roi(:,k),:),0,1);

        case 'raw'
            v(:,k) = x(roi(:,k),:);
            d = [];

        otherwise
            error('unknown option');
    end
end

if nargin>3
    v0 = image_roi(x0, roi, typ);
    v  = v - v0;
end