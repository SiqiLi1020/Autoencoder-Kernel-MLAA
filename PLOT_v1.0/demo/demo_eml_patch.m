% this is a demo for penalized likelihood PET reconstruction using 
% patch-based regularization.
%
% Guobao Wang @ UC Davis (05-07-2013)
%
clear; clc;

% remember to run "setup.m" to add the path of PLOT. If you have done so,
% please comment the following two lines.
PLOTPATH = '../';   % '../' effective only when you are in the "demo" foler 
run([PLOTPATH,'setup']);  

% image size
imgsiz = [256 256]; 

%% system matrix
sysflag = 1;  % 1: using Fessler's IRT; 
              % 0: using your own system matrix G
disp(':: 1. Generating system matrix G ...')
if sysflag
    % require Fessler's IRT matlab toolbox to generate a system matrix
    irtpath = '../../irt/setup'; % This is where your irt/setup.m sits
    run(irtpath);  

    ig = image_geom('nx', imgsiz(1), 'ny', imgsiz(2), 'fov', 32);
    
    % field of view
    ig.mask = ig.circ(10, 10) > 0;

    % system matrix G
    prjsiz = [256 192];
    sg = sino_geom('par', 'nb', 256, 'na', 192, 'dr', 32 / 256);
    G  = Gtomo2_strip(sg, ig, 'single', 1);
    Gopt.mtype = 'fessler';
    Gopt.ig    = ig;
    Gopt.imgsiz= imgsiz;
else
    Gfold = ''; % where you store your system matrix G;
    load([Gfold, 'G']);	% load matrix G
    Gopt.mtype  = 'matlab';
    Gopt.imgsiz = imgsiz;
end

%% simulate PET data
disp(':: 2. Simulating PET data ...')

% load the original image
load('brain_data.mat');

% noisy measurements
count = 2e5; % a total of 200k events
proj  = proj_forw(G, Gopt, x0); 
ai = exp(- proj_forw(G, Gopt, ct) );  % attenuation factor
ri = mean(ai(:).*proj(:)) * 0.2 * ones(size(proj));  % 20% uniform background
y0 = ai.*proj + ri; % noiseless projection
cs = count / sum(y0(:));
y0 = cs * y0;
ri = cs * ri;
yi = poissrnd(y0); % noisy projection
ni = ones(size(yi))*cs; % (uniform) normalization factor

%% reconstruction using DEM
disp(':: 3. Image reconstruction ...')

% reconstruction parameters
beta  = 2^-6;
delta = 1e-9;
maxit = 100; 
xinit = [];  % uniform initial estimate

% construct neighborhood
R = buildNbhd(imgsiz, 'clique', 2); % second order neighborhood (8 pixels)

% anisotropic regularization forms
R1 = buildPenalty(imgsiz, 'lange', delta, 1, 0, R); % pixel-based 
R2 = buildPenalty(imgsiz, 'lange', delta, 3, 0, R); % patch-based (patch size: 3x3)

% pixel-based regularization
disp('--using pixel-based regularization...')
tic; [x{1}, out{1}, Phi{1}] = eml_dem(yi, ni.*ai, G, Gopt, xinit, ri, maxit, beta, R1); toc;

% patch-based regularization 
disp('--using patch-based regularization...')
tic; [x{2}, out{2}, Phi{2}] = eml_dem(yi, ni.*ai, G, Gopt, xinit, ri, maxit, beta, R2); toc;

%% display
lname = {'Pixel, DEM', 'Patch, DEM'};
lmark = {'rx--','bo-'};

% SNR
for i = 1:length(x)
    for j = 1:size(out{i}.xest,2)
        snr{i}(j) = image_snr(x0, out{i}.xest(:,j)); 
    end
end

% reconstructed images
for i = 1:length(x)
    figure, imagesc(reshape(x{i},imgsiz), [0 100]); colormap(jet(256)); 
    axis image; colorbar; set(gca,'FontSize',16); 
    title(sprintf('%s, SNR=%3.2fdB', lname{i}, snr{i}(end)))
end

% plot SNR
figure, 
for i = 1:length(x)
    plot(out{i}.step, snr{i},lmark{i}, 'LineWidth', 3); hold on;
end
set(gca,'FontSize',16); xlabel('iteration number'); ylabel('image SNR (dB)');
legend(lname{1}, lname{2}, 'Location', 'Best')
