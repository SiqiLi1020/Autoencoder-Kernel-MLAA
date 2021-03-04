function y = proj_forw(G, Gopt, x, numfrm)
%--------------------------------------------------------------------------
% forward projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if nargin<4 | isempty(numfrm)
    numfrm = 1;
end

% frame-by-frame
if numfrm>1
    x = reshape(x,[length(x(:))/numfrm,numfrm]);
    for m = 1:numfrm
        y(:,m) = proj_forw(G, Gopt, x(:,m), 1);
    end
    return;
end

if numfrm==1
    x = x(:);
end
if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end

% use kernel
if ~isfield(Gopt,'kernel')
    Gopt.kernel = [];
end
if not(isempty(Gopt.kernel))
    x = Gopt.kernel * x;
end

switch Gopt.mtype
    case 'matlab'
        y = G * x;
        
    case 'handle'
        y = Gopt.forw(x);
        
    case 'fessler'
        y = G * x(Gopt.ig.mask,:);
        y = double(y);
        
    case 'fessler2.5'
        N = prod(Gopt.ig.dim);
        x = reshape(x,[N,Gopt.Nz]);
        y = G * x(Gopt.ig.mask,:);
        y = double(y(:));
        
    case 'usc'
        put_data('temp.image.usc', x, 'float32');
        setenv('NUMOFWORKERS', '32');
        unix(sprintf('%s %s temp.image.usc', Gopt.pforw, Gopt.cfile));
        y = get_data('temp.image.usc.forw');
        
    case 'zhou'
        put_data('temp.image', x, 'float32');
        setenv('OMP_NUM_THREADS', '32');
        unix(sprintf('%s -c %s -f temp.image -o temp.image.forw', Gopt.pforw, Gopt.cfile));
        y = get_data('temp.image.forw');
       
    case 'dst2d'
        
        acqParams.nU                   = 249;
        acqParams.nPhi                 = 210;
        acqParams.nZ                   = 47;
        acqParams.sU                   = 3.1966;
        acqParams.sV                   = 3.2646;
        scanner.numBlocksPerRing       = 70;
        scanner.radialCrystalsPerBlock = 6;
        scanner.radBlockSize           = 38.35;
        subsetList                     = 0;
        reconParams.FOV                = 500;
        reconParams.nx                 = 128;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = -0.17089;
        reconParams.yOffset            = 1.6078;
        reconParams.rotate             = -4.6979;
        scanner.effectiveRingDiameter  = 903;
        
        x = reshape(x,reconParams.nx, reconParams.nx, acqParams.nZ);
        y = FDD2D_ng(single(x), acqParams.nU, acqParams.nPhi, reconParams.nx, ...
                     acqParams.nZ, acqParams.sU, reconParams.FOV/reconParams.nx, ...
                     acqParams.sV, scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, ...
                     scanner.radBlockSize, subsetList, reconParams.numSubsets, ...
                     reconParams.xOffset, reconParams.yOffset, reconParams.rotate, ...
                     scanner.effectiveRingDiameter);
        y = double(y(:));  
        
    case 'gedst'
        
        acqParams.nU                   = 249;
        acqParams.nV                   = 553;
        acqParams.nPhi                 = 210;
        acqParams.nZ                   = 47;
        acqParams.sU                   = 3.1966;
        acqParams.sV                   = 3.2646;
        scanner.numBlocksPerRing       = 70;
        scanner.radialCrystalsPerBlock = 6;
        scanner.radBlockSize           = 38.35;
        subsetList                     = 0;
        reconParams.FOV                = 500;
        reconParams.nx                 = 128;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = -0.17089;
        reconParams.yOffset            = 1.6078;
        reconParams.rotate             = -4.6979;
        scanner.effectiveRingDiameter  = 903;
        
        x = reshape(x,reconParams.nx, reconParams.nx, acqParams.nZ);
        y = FDD3D_ng(single(x), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        y = double(y(:));
       
    case '6902d'
        
        acqParams.nU                   = 381;
        acqParams.nPhi                 = 288;
        acqParams.nZ                   = 1;
        acqParams.sU                   = 2.1306;
        acqParams.sV                   = 3.27;
        scanner.numBlocksPerRing       = 64;
        scanner.radialCrystalsPerBlock = 9;
        scanner.radBlockSize           = 38.35;
        subsetList                     = 0;
        reconParams.FOV                = 700;
        reconParams.nx                 = 192;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = 0;
        reconParams.yOffset            = 0;
        reconParams.rotate             = -5.02;
        scanner.effectiveRingDiameter  = 829;
        
        x = reshape(x,reconParams.nx, reconParams.nx, acqParams.nZ);
        y = FDD2D_ng(single(x), acqParams.nU, acqParams.nPhi, reconParams.nx, ...
                     acqParams.nZ, acqParams.sU, reconParams.FOV/reconParams.nx, ...
                     acqParams.sV, scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, ...
                     scanner.radBlockSize, subsetList, reconParams.numSubsets, ...
                     reconParams.xOffset, reconParams.yOffset, reconParams.rotate, ...
                     scanner.effectiveRingDiameter);
        y = double(y(:));  
        
    case 'ge690'

        acqParams.nU                   = 381;
        acqParams.nV                   = 553;
        acqParams.nPhi                 = 288;
        acqParams.nZ                   = 47;
        acqParams.sU                   = 2.1306;
        acqParams.sV                   = 3.27;
        scanner.numBlocksPerRing       = 64;
        scanner.radialCrystalsPerBlock = 9;
        scanner.radBlockSize           = 38.35;
        subsetList                     = 0;
        reconParams.FOV                = 700;
        reconParams.nx                 = 192;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = 0;
        reconParams.yOffset            = 0;
        reconParams.rotate             = -5.02;
        scanner.effectiveRingDiameter  = 829;
        
        x = reshape(x,reconParams.nx, reconParams.nx, acqParams.nZ);
        y = FDD3D_ng(single(x), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        y = double(y(:));  
        
end
