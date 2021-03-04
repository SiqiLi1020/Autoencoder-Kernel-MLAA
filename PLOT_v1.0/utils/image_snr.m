function s = snr(x0, x)

s = -10*log10( sum((x0(:)-x(:)).^2)/sum(x0(:).^2) );