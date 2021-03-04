function mip_orthview(x, maxval, color_scale, pos, window, flag, printname)
% display 3D images
%
% gbwang@ucdavis.edu (01-09-2013)
%
if flag==1
    pos(2) = size(x,2)-pos(2);
end
if nargin<2 | isempty(maxval)
    maxval = max(x(:));
end
if nargin<3 | isempty(color_scale)
    color_scale = 'jet';
end
thresh = [0 maxval];
if nargin<5 | isempty(window)
    window = [1 size(x,1) 1 size(x,2) 1 size(x,3)];
end
if nargin<6
    flag = 0;
end
if nargin<7
    printname = [];
end
if window(6)>size(x,3)
    x(:,:,(size(x,3)+1):window(6)) = 0;
end
if nargin<4 | isempty(pos)
    x = x(window(1):window(2), window(3):window(4), window(5):window(6));
    x1 = squeeze(max(x,[],1));
    x2 = squeeze(max(x,[],2));
    x3 = squeeze(max(x,[],3));
    orth = 0;
else
    x1 = squeeze(x(pos(1),:,:));
    x2 = squeeze(x(:,pos(2),:));
    x3 = squeeze(x(:,:,pos(3)));
    orth = 1;
    x1 = x1(window(3):window(4), window(5):window(6));
    x2 = x2(window(1):window(2), window(5):window(6));
    x3 = x3(window(1):window(2), window(3):window(4));
    x = x(window(1):window(2), window(3):window(4), window(5):window(6));
end
    
figure,
if flag
    imagesc(flipud(imrotate(x1,-90)), thresh); 
else
    imagesc(flipud(imrotate(x1,-90)), thresh); 
end
colormap(color_scale);
axis image; 
axis off;
%title('Coronal');
if orth
    hold on; 
    p = window(4)-pos(2);
    line([p,p],[1 5],'color','r');
    line([p,p],[size(x1,2)-5 size(x1,2)],'color','r');
    p = -pos(3)+window(6)+1;
    line([1 5],[p, p],'color','r');
    line([size(x1,1)-5 size(x1,1)],[p, p],'color','r');
end
if not(isempty(printname))& ~flag
    colorbar;set(gca, 'FontSize', 16);
    print([printname,'_sag.png'],'-dpng')
end

figure,
if flag
    imagesc(imrotate(x2,90), thresh);
else
    imagesc(fliplr(imrotate(x2,90)), thresh);
end
colormap(color_scale);
axis image; 
axis off;
%title('Sagital')
if orth
    hold on; 
    p = size(x2,1)-(pos(1)-window(1))+1;
    line([p,p],[1 5],'color','r');
    line([p,p],[size(x2,2)-5 size(x2,2)],'color','r');
    p = -pos(3)+window(6)+1;
    line([1 5],[p, p],'color','r');
    line([size(x2,1)-5 size(x2,1)],[p, p],'color','r');
end
if not(isempty(printname))& ~flag
    colorbar; set(gca, 'FontSize', 16);
    print([printname,'_cor.png'],'-dpng')
end

figure,
if flag
    imagesc(imrotate(x3,-90), thresh); 
else
    imagesc(flipud(x3), thresh); 
end
axis image;
colormap(color_scale);
axis image; 
axis off;
axis image;
%title('Transverse');
if orth
    hold on; 
    p = window(2)-pos(1)+1;
    line([1 5],[p, p],'color','r');
    line([size(x3,2)-5 size(x3,2)],[p, p],'color','r');
    p = pos(2)-window(3)+1;
    line([p, p], [1 5],'color','r');
    line([p, p], [size(x3,1)-5 size(x3,1)],'color','r');
end
if not(isempty(printname))
    colorbar; set(gca, 'FontSize', 16);
    print([printname,'_tra.png'],'-dpng')
end

