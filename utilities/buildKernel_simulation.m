function [N, W] = buildKernel_simulation(imgsiz, nbrtyp, nbrpar, X, kertyp, kerpar, normflag, midflag, inplaneflag)
% gbwang@ucdavis.edu (11-01-2012)
imgdim = length(find(imgsiz>1));
if imgdim==1
    imgsiz = [imgsiz(imgsiz(:)>1) 1 1];
elseif imgdim==2
    imgsiz = [imgsiz(imgsiz(:)>1) 1];
end
if nargin<2 | isempty(nbrtyp)
    nbrtyp = 'cube';
end
if nargin<3 | isempty(nbrpar)
    nbrpar = 3;
end
if length(nbrpar)==1
    if imgdim==3
        nbrpar = [nbrpar nbrpar nbrpar];
    elseif imgdim==2
        nbrpar = [nbrpar nbrpar 1];
    end
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
if nargin<8 | isempty(midflag)
    midflag = 1;
end
if nargin<9 | isempty(inplaneflag)
    inplaneflag = 0;
end
numvox = prod(imgsiz);
J = [1:numvox]';
switch nbrtyp
    case 'clique'
        s   = [1  0  1  1  0  1  0 -1  0  1 -1 -1  1;
            0  1  1 -1  0  0  1  0 -1  1  1 -1 -1;
            0  0  0  0  1  1  1  1  1  1  1  1  1];
        d   = sqrt(sum(abs(s).^2,1));
        if imgdim==1
            idx = 1;
        elseif imgdim==2
            idx = 1:4;
        elseif imgdim==3
            idx = 1:13;
        end
        s = s(:,idx); d = d(idx);
        if nbrpar==0.5
            idx = d==1;
        elseif nbrpar==1
            idx = d==1;
        elseif nbrpar==2
            idx = d>0;
        end
        s = s(:,idx); d = d(idx);
        N = zeros(numvox, size(s,2));
        W = zeros(size(N));
        for k = 1:size(s,2)
            N(:,k) = setBoundary3_simulation(J, s(:,k), imgsiz);
            W(:,k) = 1/d(k);
        end
        if nbrpar>=1
            for k = 1:size(s,2)
                N(:,k+size(s,2)) = setBoundary3_simulation(J, -s(:,k), imgsiz);
                W(:,k+size(s,2)) = 1/d(k);
            end
        end
        if midflag
            N = [J N];
            W = [sum(W,2)/4 W];
        end
        if not(isempty(X))
            W = W.*calc_wgt_simulation(X, N, kertyp, kerpar);
        end
    case 'cube'
        wlen = 2*floor(nbrpar/2);
        xidx = -wlen(1)/2:wlen(1)/2;
        yidx = -wlen(2)/2:wlen(2)/2;
        zidx = -wlen(3)/2:wlen(3)/2;
        if imgsiz(3)==1 zidx = 0; end
        [I1, I2, I3] = ndgrid(1:imgsiz(1),1:imgsiz(2),1:imgsiz(3));
        if imgdim==2
            h = fspecial('gaussian', nbrpar(1), nbrpar(1)/(4*sqrt(2*log(2))));
        elseif imgdim==3
            h = fspecial3_simulation('gaussian', nbrpar);
        end
        N = zeros(numvox, length(h(:)));
        W = N;
        l = 1; n = 1;
        for x = xidx
            Xnew = setBoundary1_simulation(I1 + x, imgsiz(1));
            for y = yidx
                Ynew = setBoundary1_simulation(I2 + y, imgsiz(2));
                for z = zidx
                    Znew = setBoundary1_simulation(I3 + z, imgsiz(3));
                    if inplaneflag & (abs(x)>0&abs(y)>0&abs(z)>0)
                        disp(sprintf('skip %d',n));
                        n = n + 1;
                        continue;
                    end
                    N(:,l) = Xnew + (Ynew-1).*imgsiz(1) + (Znew-1)*imgsiz(1)*imgsiz(2);
                    W(:,l) = h(n);
                    l = l + 1;
                    n = n + 1;
                end
            end
        end
        if l<=size(N,2)
            N(:,l:end) = [];
            W(:,l:end) = [];
        end
        if ~midflag
            i = ceil(size(N,2)/2);
            N(:,i) = [];
            W(:,i) = [];
        end
        if not(isempty(X))
            W = W.*calc_wgt_simulation(X, N, kertyp, kerpar);
        end
    case 'knn'
        k = nbrpar(1);
        [N, D] = knnsearch(X, X, 'dist', 'seuclidean', 'k', k);
        if ~midflag
            N = N(:,2:end);
        end
        W = calc_wgt_simulation(X, N, kertyp, kerpar);
end
if normflag
    sumw = repmat(sum(W,2),[1 size(W,2)]);
    W = W./sumw;
    W(sumw==0) = 0;
end

function J = setBoundary3_simulation(J, d, imgsiz)
[x,y,z] = ind2sub(imgsiz,J);
x = setBoundary1_simulation(x+d(1), imgsiz(1));
y = setBoundary1_simulation(y+d(2), imgsiz(2));
z = setBoundary1_simulation(z+d(3), imgsiz(3));
J = sub2ind(imgsiz,x,y,z);
function x = setBoundary1_simulation(x, N)
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
function W = calc_wgt_simulation(X, N, kertyp, kerpar)
J  = [1:size(N,1)]';
D  = zeros(size(N));
for i = 1:size(N,2)
    D(:,i) = sqrt(mean(((X(J,:)-X(N(:,i),:))*sparse(1:size(X,2),1:size(X,2),1./std(X,0,1))).^2,2));
end
switch kertyp
    case 'invdist'
        W = 1./D;
    case 'radial'
        W = exp(-D.^2./(2*kerpar.^2));
    case 'poly'
        for i = 1:size(N,2)
            W(:,i) = ( mean(X(J,:).*X(N(:,i),:),2) + kerpar ).^2;
        end
    case 'nnlle'
        for j = 1:size(N,1)
            w = lsqnonneg(X(N(j,:),:)',X(j,:)');
            if any(isnan(w))
                w = zeros(length(w),1);
            end
            W(j,:) = w;
        end
    otherwise
        error('unknown kernel type')
end
function h = fspecial3_simulation(type,siz)
error(nargchk(1,2,nargin))
type = lower(type);

if nargin==1
    siz = 5;
end
if numel(siz)==1
    siz = round(repmat(siz,1,3));
elseif numel(siz)~=3
    error('Number of elements in SIZ must be 1 or 3')
else
    siz = round(siz(:)');
end
switch type
    case 'average'
        h = ones(siz)/prod(siz);
    case 'gaussian'
        sig = siz/(4*sqrt(2*log(2)));
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
    case 'ellipsoid'
        R = siz/2;
        R(R==0) = 1;
        h = ones(siz);
        siz = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        I = (x.*x/R(1)^2+y.*y/R(2)^2+z.*z/R(3)^2)>1;
        h(I) = 0;
        h = h/sum(h(:));
    case 'laplacian'
        h = zeros(3,3,3);
        h(:,:,1) = [0 3 0;3 10 3;0 3 0];
        h(:,:,3) = h(:,:,1);
        h(:,:,2) = [3 10 3;10 -96 10;3 10 3];
    case 'log'
        sig = siz/(4*sqrt(2*log(2)));
        siz   = (siz-1)/2;
        [x,y,z] = ndgrid(-siz(1):siz(1),-siz(2):siz(2),-siz(3):siz(3));
        h = exp(-(x.*x/2/sig(1)^2 + y.*y/2/sig(2)^2 + z.*z/2/sig(3)^2));
        h = h/sum(h(:));
        arg = (x.*x/sig(1)^4 + y.*y/sig(2)^4 + z.*z/sig(3)^4 - ...
            (1/sig(1)^2 + 1/sig(2)^2 + 1/sig(3)^2));
        h = arg.*h;
        h = h-sum(h(:))/prod(2*siz+1);
    otherwise
        error('Unknown filter type.')
end