function [v, d] = image_crc(x, roi, typ, bkg, bkgID)
% tumor contrast recovery
%
% gbwang@ucdavis.edu Jan. 2012
%

if isvector(x)
    x = x(:);
end
if nargin<5 | isempty(bkgID)
    temp = unique(bkg(:));
    bkgID = temp(temp>10);
end

roi = roi>0;
switch typ
    case 'crc'
        roimean = mean(x(roi,:));
        bkgmean = mean(x(bkg,:)); 
        v = (roimean-bkgmean)./bkgmean;
        d = std(x(bkg,:),0,1)./bkgmean*100;
        
    case 'nema'
        roimean = mean(x(roi,:));
        for i = 1:length(bkgID)
            xi(i,:) = mean(x(bkg==bkgID(i),:));
        end
        bkgmean = mean(xi,1);
        v = (roimean-bkgmean)./bkgmean;
        d = std(xi,0,1)./bkgmean*100;
        
    otherwise
        error('unknown option');
    end

end