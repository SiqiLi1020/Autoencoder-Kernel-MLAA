% A testing demo: Modified Kernel MLAA using Autoencoder for PET-enabled Dual-Energy CT. 
% The details  are described in
% Siqi Li and Guobao Wang. "Modified Kernel MLAA using Autoencoder for PET-enabled Dual-Energy CT." arXiv preprint arXiv: https://arxiv.org/abs/2010.07484 (2020).
%-----------------------------------------------------------------------------------------------------------------------------------------------------------------
% Programmer: Siqi Li and Guobao Wang, UC DAVIS.
% Contact: sqlli@ucdavis.edu


clc;
clear;
i = 1;

% count level
count = 5e6;

% choose rectype 'mlaa','standard kernel mlaa','RED kernel mlaa','Unet kernel mlaa';
rectype = 'Unet kernel mlaa';

% initialization method
initype = 'CT';

run('../KER_v0.11/KER_v0.11/setup');  
run('../PLOT_v1.2/setup');
addpath('utils');

%% load data

% load the image data at 511 keV
load('data/psct_xcat');

% load the correction data
load(sprintf('data/proj%d',0));

% system matrix - attenuation
load('data/GE690TOF2D/GE690_attn_sysmat_180x3.27mm.mat');
Gopt.imgsiz = [180 180];
Gopt.prjsiz = [281 288];

% system matrix - emission
load('data/GE690TOF2D/GE690_tof2d_sysmat_11tbin_180x3.27mm.mat');
Popt.imgsiz = [180 180];
numbin      = 11;
Popt.prjsiz = [281 288 numbin];

% joint
Gopt.attnFlag = 1;
Popt.emisFlag = 1;
Popt.disp     = 1;
Popt.savestep = 100;
maxit = 3000;

% mask 
Gopt.mask = mask; 

% GCT initial
imgsiz = size(mask);
uinit = CT2LAC(CT,'120','bilinear')/10;
inilabel = '';
if strcmp(initype,'uniform')
    uinit = zeros(imgsiz);
    uinit(mask) = 0.01;
    inilabel = '-UI';
end

% PET activity initial
xinit = zeros(size(mask));
xinit(Gopt.mask) = 1;

% fold
mkdir('result/', sprintf('xcat_proj%dm_rec',count/1e6));
%foldx = sprintf('../data/xcat_proj%dm_rec',count/1e6);

% load noisy projection
load(sprintf('data/proj%d', i));

% %PET initial
if ~strcmp(rectype,'mlem')
    load(sprintf('data/%s_%d','mlem', i));
    xinit = out.xest(:,end);
end


%% different kernel matrix

% standard kernel matrix
imgsiz_CT = size(CT);
R = buildNbhd(imgsiz_CT, 'clique', 1); % Extract features using a 3x3 patch 
I = [[1:prod(imgsiz_CT)]' R.N]; 
F = CT(I);
F = F * diag(1./std(F,1)); % normalization
sigma = 1;
[N, W] = buildKernel(imgsiz, 'knn', 50, F, 'radial', sigma, 1);
K = buildSparseK(N, W);

% RED-CNN kernel matrix
load 'trained feature maps/red_fe';
F_fe = permute(red1_fe_96,[2,3,1]); 
f_map1 = zeros(180*180,size(F_fe,3));
for m = 1: size(F_fe,3)
     feature_map = imrotate(flipud(F_fe(:,:,m)), -90);
     f_map1(:,m) = feature_map(:);
end
f_map1 = f_map1 * diag(1./std(f_map1,1));
sigma = 1;
[N, W] = buildKernel(imgsiz, 'knn', 50, f_map1, 'radial', sigma, 1);
K_RED = buildSparseK(N, W);

% Unet kernel matrix
load 'trained feature maps/U_fe';
F_fe = permute(U_fe3,[2,3,1]);
f_map1 = zeros(180*180,size(F_fe,3));
 for m = 1 : size(F_fe,3)
    feature_map = imrotate(flipud(F_fe(:,:,m)), -90);
    f_map1(:,m) = feature_map(:);
end
f_map1 = f_map1 * diag(1./std(f_map1,1));
sigma = 1;
[N, W] = buildKernel(imgsiz, 'knn', 50, f_map1, 'radial', sigma, 1);
K_Unet = buildSparseK(N, W);

%% Synergistic reconstruction
switch rectype
    
    case 'mlem'
        Popt.savestep = 10;
        a0 = exp(- proj_forw(G, Gopt, uinit) );  % attenuation factor
        ai = repmat(a0,[1 Popt.prjsiz(3)]);
        [x, out] = eml_dem(yi, ni.*ai(:), P, Popt, xinit(:), ri, 100);
    
    case 'mlaa'
        [u, x, out] = psct_kmlaa(yi, ni, G, Gopt, uinit(:), P, Popt, xinit(:), ri, maxit);
        save(sprintf('result/xcat_proj%dm_rec/%s_%d',count/1e6,'mlaa', i),'out');
        
    case 'standard kernel mlaa'
        [u, x, out] = psct_kmlaa(yi, ni, G, Gopt, uinit(:), P, Popt, xinit(:), ri, maxit, K);
        out.uest = K * out.uest;
        save(sprintf('result/xcat_proj%dm_rec/%s_%d',count/1e6,'kmlaa', i),'out');

    case 'RED kernel mlaa'
        [u, x, out] = psct_kmlaa(yi, ni, G, Gopt, uinit(:), P, Popt, xinit(:), ri, maxit, K_RED);
        out.uest = K_RED * out.uest;
        save(sprintf('result/xcat_proj%dm_rec/%s_%d',count/1e6,'kmlaa-red', i),'out');
        
    case 'Unet kernel mlaa'
        [u, x, out] = psct_kmlaa(yi, ni, G, Gopt, uinit(:), P, Popt, xinit(:), ri, maxit, K_Unet);
        out.uest = K_Unet * out.uest;
        save(sprintf('result/xcat_proj%dm_rec/%s_%d',count/1e6,'kmlaa-Unet', i),'out');

end


