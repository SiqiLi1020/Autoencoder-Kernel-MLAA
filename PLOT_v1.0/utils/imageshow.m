function imageshow(x, imgsiz, clrbar, ang, clrmap)

if nargin<2 | isempty(imgsiz)
    imgsiz = size(x);
end
if nargin<3 | isempty(clrbar)
    clrbar = [min(x(:)), max(x(:))];
end
if nargin<4 | isempty(ang)
    ang = -90;
end
if nargin<5 | isempty(clrmap)
    clrmap = 'jet';
end

imagesc(fliplr(imrotate(reshape(x,imgsiz),ang)), clrbar); 
axis image; 
colorbar;