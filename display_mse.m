clc;
clear;

load ('data/psct_xcat');

load 'result/xcat_proj5m_rec/kmlaa_1';
KMLAA = out;

load 'result/xcat_proj5m_rec/kmlaa-Unet_1';
DIP = out;

load 'result/xcat_proj5m_rec/mlaa_1';
Mlaa = out;

load 'result/xcat_proj5m_rec/kmlaa-red_1';
Red = out;

for i = 1:31
    Redmse(i) = 10*log10(sum((Red.uest(:,i).*10 - u0(:)).^2) / sum(u0(:).^2));
    KMLAAmse(i) = 10*log10(sum((KMLAA.uest(:,i).*10 - u0(:)).^2) / sum(u0(:).^2));
    DIPmse(i) = 10*log10(sum((DIP.uest(:,i).*10 - u0(:)).^2) / sum(u0(:).^2));
    MLAAmse(i) = 10*log10(sum((Mlaa.uest(:,i).*10 - u0(:)).^2) / sum(u0(:).^2));
end

x = 0:100:3000;
x(1) = 1;
y=1:31;
plot(x, MLAAmse(y),'m*-',x,KMLAAmse(y), 'rx-', x,Redmse(y),'gs-',x,DIPmse(y),'bo-','LineWidth', 2);
ylim([-29 -8]);
xlabel('iteration number', 'FontName','Times New Roman', 'FontWeight','normal','Fontsize',13);
ylabel('MSE ','FontName','Times New Roman', 'FontWeight','normal','Fontsize',13);
%title(sprintf('MSE'),'FontName','Times New Roman', 'FontWeight','normal','Fontsize',13); 
h = legend('Standard MLAA','Standard kernel MLAA', 'RED kernel MLAA','Unet kernel MLAA','Location','NorthEast');
set(h,'FontName','Times New Roman', 'FontWeight','normal','Fontsize',13);