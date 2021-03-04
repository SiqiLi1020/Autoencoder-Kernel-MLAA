function H = buildNbhd(imgsiz, nbrtyp, nbrpar, X, kertyp, kerpar, normflag)
%
% Get the index and weight of neighboring voxes in a neighborhood
%
% Guobao Wang @ UC Davis (10-01-2012)
%
% last modified: 05-20-2013
%

%% check input variables
imgdim = length(find(imgsiz>1));
if imgdim==1
    imgsiz = [imgsiz(imgsiz(:)>1) 1 1];
elseif imgdim==2
    imgsiz = [imgsiz(imgsiz(:)>1) 1];
end
if nargin<2 | isempty(nbrtyp)
    nbrtyp = 'clique';
end
if nargin<3 | isempty(nbrpar)
    nbrpar = 2;
end
if nargin<4
    X = [];
end
if nargin<5 | isempty(kertyp)
    kertyp = 'invdist';
end
if nargin<6 | isempty(kerpar)
    kerpar = 1;
end
if nargin<7 | isempty(normflag)
    normflag = 1;
end

numvox = prod(imgsiz);
J = [1:numvox]';

%% neighborhood type
switch nbrtyp
    
    case 'clique'
       
        s   = [1  0  1  1  0  1  0 -1  0  1 -1 -1  1;
               0  1  1 -1  0  0  1  0 -1  1  1 -1 -1;
               0  0  0  0  1  1  1  1  1  1  1  1  1];
        d   = sqrt(sum(abs(s).^2,1));

        % image dimension
        if imgdim==1
            idx = 1;
        elseif imgdim==2
            idx = 1:4;
        elseif imgdim==3
            idx = 1:13;
        end
        s = s(:,idx); d = d(idx);
        
        % neighborhood order
        if nbrpar==0.5  % for TV
            idx = d==1;
        elseif nbrpar==1
            idx = d==1;
        elseif nbrpar==2
            idx = d>0;
        end
        s = s(:,idx); d = d(idx);
        
        % neighboring voxel index
        N = zeros(numvox, size(s,2));
        D = zeros(size(N));
        for k = 1:size(s,2)
            N(:,k) = setBoundary3(J, s(:,k), imgsiz);
            D(:,k) = d(k);
        end
        if nbrpar>=1
            for k = 1:size(s,2)
                N(:,k+size(s,2)) = setBoundary3(J, -s(:,k), imgsiz);
                D(:,k+size(s,2)) = d(k);
            end
        end
        
    case 'cube'
        
        % sizes of neighborhood window
        wlen = 2*floor(nbrpar/2); % length of neighborhood window
        widx = -wlen/2:wlen/2; 
        xidx = widx; yidx = widx; zidx = widx;
        if imgsiz(3)==1 zidx = 0; end
        
        % image grid
        [I1, I2, I3] = ndgrid(1:imgsiz(1),1:imgsiz(2),1:imgsiz(3));
        
        % index and distance
        N = zeros(numvox, length(widx));
        D = zeros(size(N));
        l = 1;
        for x = xidx
            Xnew = setBoundary1(I1 + x, imgsiz(1));
            for y = yidx
                Ynew = setBoundary1(I2 + y, imgsiz(2));
                for z = zidx
                    Znew = setBoundary1(I3 + z, imgsiz(3));
                    N(:,l) = Xnew + (Ynew-1).*imgsiz(1) + (Znew-1)*imgsiz(1)*imgsiz(2);
                    D(:,l) = sqrt(x^2+y^2+z^2);
                    l = l + 1;
                end
            end
        end
        
        % removing the middle column
        i = floor((size(N,2)+1)/2);
        N(:,i) = []; 
        D(:,i) = [];

    case 'knn'
        
        k = nbrpar;
        N = knnsearch(X, X, 'dist', 'seuclidean', 'k', k+1);
        N = N(:,2:end);

end

% distance
if not(isempty(X))
    sd = diag(1./std(X,0,1));
    D = zeros(size(N));
    for i = 1:size(N,2)
        D(:,i) = sqrt(sum(((X(J,:)-X(N(:,i),:))*sd).^2,2));
    end
end

% kernel type
switch kertyp
    case 'invdist'
        W = 1./D;
        
    case 'radial' % radial Gaussian
        W = exp(-D.^2/kerpar);
        
    case 'poly' % polynoimial
        for i = 1:size(N,2)
            W(:,i) = ( sum(X(J,:).*X(N(:,i),:),2) + kerpar ).^2;
        end
        
    otherwise
        error('unknown kernel type')
end

% normalized to 1 and output
H.imgsiz = imgsiz;
if normflag
    H.W = W./repmat(sum(W,2),[1 size(W,2)]);
else
    H.W = W;
end
H.N = N;
H.nbrtyp = nbrtyp;
H.nbrpar = nbrpar;
H.kertyp = kertyp;
H.kerpar = kerpar;

%% create a sparse matrix C 
H.C = buildSparseCK(H.N, H.W, J, numvox);


%% sub functions

%--------------------------------------------------------------------------
function J = setBoundary3(J, d, imgsiz)
%--------------------------------------------------------------------------
[x,y,z] = ind2sub(imgsiz,J);
x = setBoundary1(x+d(1), imgsiz(1)); 
y = setBoundary1(y+d(2), imgsiz(2)); 
z = setBoundary1(z+d(3), imgsiz(3));
J = sub2ind(imgsiz,x,y,z);
    
%--------------------------------------------------------------------------
function x = setBoundary1(x, N)
%--------------------------------------------------------------------------
if N==1
    x = x(:).^0;
else
    idx = x(:)>N;
    if any(idx)
        x(idx) = N - (x(idx)-N); 
    end
    idx = x(:)<1;
    x(idx) = 1 + (1-x(idx)); 
    x = x(:);
end

