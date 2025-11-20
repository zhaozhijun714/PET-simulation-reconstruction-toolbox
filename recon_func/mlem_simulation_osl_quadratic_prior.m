function x = mlem_simulation_osl_quadratic_prior(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, beta, Ndx, Ndy)
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end
if isempty(maxit)
    maxit = 10;
end
sens = proj_back_simulation(G, ni, imgsiz, voxel_size, use_psf, FWHM);
mask = sens > 0;
x    = max(mean(x0(:)) * 1e-9, x0(:));
x(~mask) = 0;
yeps = mean(yi(:)) * 1e-9;
wx   = sens;
weights = compute_weights_simulation(voxel_size, Ndx, Ndy);
weights_quad = quadWeights_simulation(weights, Ndx, Ndy);
for it = 1 : maxit
    disp(sprintf('iteration %d',it));
    grad = Quadratic_prior_simulation(x, weights_quad, imgsiz(1), imgsiz(2), Ndx, Ndy);
    wx_osl = wx + beta * grad;
    wx_osl(wx_osl < 0) = 0;
    yb = ni .* proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM) + ri;
    yy = yi ./ (yb + yeps);
    yy(yb == 0 & yi == 0) = 1;
    xb = proj_back_simulation(G, ni.*yy, imgsiz, voxel_size, use_psf, FWHM);
    x  = x ./ wx_osl .* xb;
    x(isinf(x)) = 0;
    x(isnan(x)) = 0;
    x(~mask) = 0;
end
end
function grad = Quadratic_prior_simulation(x, weights_quad, Nx, Ny, Ndx, Ndy)
x = x(:);
x = padding_simulation(reshape(x, Nx ,Ny), [Ndx Ndy]);
weights_quad = -weights_quad;
weights_quad(ceil(numel(weights_quad) / 2)) = -weights_quad(ceil(numel(weights_quad) / 2));
grad = convn(x, weights_quad, 'valid');
grad(isinf(grad)) = 0;
grad(isnan(grad)) = 0;
grad = grad(:);
end
function weights_quad = quadWeights_simulation(weights, Ndx, Ndy)
weights_quad = weights / sum(weights(~isinf(weights)));
weights_quad = [weights_quad(1 : floor(end / 2)); abs(sum(weights_quad(~isinf(weights)))); weights_quad(ceil(end / 2) + 1 : end)];
weights_quad = reshape(weights_quad, Ndx * 2 + 1, Ndy * 2 + 1);
end
function A = padding_simulation(A,sizeP,varargin)
if nargin == 2 || (nargin >= 3 && (isempty(varargin{1}) || ~strcmp(varargin{1},'zeros')))
    A = [flipud(A(1:sizeP(2),:,:,:));A;flipud(A(end-sizeP(2) + 1:end,:,:,:))];
    A = [fliplr(A(:,1:sizeP(1),:,:)),A,fliplr(A(:,end-sizeP(1) + 1:end,:,:))];
    if length(sizeP) == 3 && sizeP(3) ~= 0
        A = cat(3, flip(A(:,:,1:sizeP(3),:),3), A);
        A = cat(3, A, flip(A(:,:,end - sizeP(3) + 1: end,:),3));
    end
else
    [x, y , z] = size(A);
    A = [zeros(sizeP(1),y + sizeP(2)*2,z);zeros(x,sizeP(2),z),A,zeros(x,sizeP(2),z);zeros(sizeP(1),y+sizeP(2)*2,z)];
    if length(sizeP) == 3 && sizeP(3) ~= 0
        A = cat(3, zeros(size(A,1),size(A,2),sizeP(3)), A);
        A = cat(3, A, zeros(size(A,1),size(A,2),sizeP(3)));
    end
end
end