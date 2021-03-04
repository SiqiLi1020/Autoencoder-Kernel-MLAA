function x = proj_back(G, Gopt, y, numfrm)
%--------------------------------------------------------------------------
% back projection 
%
% Guobao Wang @ UC Davis (10-01-2012)
%

if nargin<4 | isempty(numfrm)
    numfrm = 1;
end

% frame-by-frame
if numfrm>1
    y = reshape(y,[length(y(:))/numfrm,numfrm]);
    
    for m = 1:numfrm
        x(:,m) = proj_back(G, Gopt, y(:,m), 1);
    end
    return;
end

if numfrm==1
    y = y(:);
end
if ~isfield(Gopt,'mtype')
    Gopt.mtype = 'matlab';
end



if ~isfield(Gopt,'kernel')
    Gopt.kernel = [];
end

switch Gopt.mtype
    case 'matlab'
        x = G' * y;
        
    case 'handle'
        x = Gopt.back(y);
        
    case 'fessler'
        a = double( G' * y ); 
        for m = 1:size(a,2)
            temp = Gopt.ig.embed(a(:,m));
            x(:,m) = temp(:);
        end
        
    case 'fessler2.5'
        M = Gopt.sg.na * Gopt.sg.nb;
        y = reshape(y,[M,Gopt.Nz]);
        a = G' * y; 
        for m = 1:size(a,2)
            temp = Gopt.ig.embed(a(:,m));
            x(:,m) = temp(:);
        end
        x = double(x(:));
        
    case 'usc'
        put_data('temp.proj.usc', y, 'float32');
        setenv('NUMOFWORKERS', '32');
        unix(sprintf('%s %s temp.proj.usc', Gopt.pback, Gopt.cfile));
        x = get_data('temp.proj.usc.back');
        
    case 'zhou'
        put_data('temp.proj', y, 'float32');
        setenv('OMP_NUM_THREADS', '32');
        unix(sprintf('%s -c %s -b temp.proj -o temp.proj.back', Gopt.pback, Gopt.cfile));
        x = get_data('temp.proj.back');
        
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
        
        y = reshape(y, acqParams.nU, acqParams.nPhi, acqParams.nZ);
        x = BDD2D_ng(single(y), acqParams.nU, acqParams.nPhi, reconParams.nx, ...
                     acqParams.nZ, acqParams.sU, reconParams.FOV/reconParams.nx, ...
                     acqParams.sV, scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, ...
                     scanner.radBlockSize, subsetList, reconParams.numSubsets, ...
                     reconParams.xOffset, reconParams.yOffset, reconParams.rotate, ...
                     scanner.effectiveRingDiameter);
        x = double(x(:));      

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
        reconParams.FOV                = 700;
        reconParams.nx                 = 128;
        reconParams.numSubsets         = 1;
        reconParams.xOffset            = -0.17089;
        reconParams.yOffset            = 1.6078;
        reconParams.rotate             = -4.6979;
        scanner.effectiveRingDiameter  = 903;
        
        y = reshape(y, acqParams.nU, acqParams.nV, acqParams.nPhi);
        x = BDD3D_ng(single(y), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        x = double(x(:));
       
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
        
        y = reshape(y, acqParams.nU, acqParams.nPhi, acqParams.nZ);
        x = BDD2D_ng(single(y), acqParams.nU, acqParams.nPhi, reconParams.nx, ...
                     acqParams.nZ, acqParams.sU, reconParams.FOV/reconParams.nx, ...
                     acqParams.sV, scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, ...
                     scanner.radBlockSize, subsetList, reconParams.numSubsets, ...
                     reconParams.xOffset, reconParams.yOffset, reconParams.rotate, ...
                     scanner.effectiveRingDiameter);
        x = double(x(:));      

        
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
        
        y = reshape(y, acqParams.nU, acqParams.nV, acqParams.nPhi);
        x = BDD3D_ng(single(y), acqParams.nU, acqParams.nV, acqParams.nPhi, acqParams.sU, acqParams.sV, ...
                     reconParams.nx, acqParams.nZ, reconParams.FOV/reconParams.nx, acqParams.sV, ...
                     scanner.numBlocksPerRing, scanner.radialCrystalsPerBlock, scanner.radBlockSize, ...
                     subsetList, reconParams.numSubsets, reconParams.xOffset, reconParams.yOffset, ...
                     reconParams.rotate, scanner.effectiveRingDiameter);
        x = double(x(:));
end   

if not(isempty(Gopt.kernel))
    x = Gopt.kernel' * x;
end
