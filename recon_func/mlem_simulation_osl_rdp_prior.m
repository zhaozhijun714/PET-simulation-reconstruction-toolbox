function x = mlem_simulation_osl_rdp_prior(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, beta, RDP_gamma, Ndx, Ndy)
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
neighbour_offsets = compute_offsets_simulation(imgsiz, Ndx, Ndy);
weights = compute_weights_simulation(voxel_size, Ndx, Ndy);
weights_RDP = RDP_weights_simulation(weights);
for it = 1 : maxit
    disp(sprintf('iteration %d',it));
    grad = RDP_simulation(x, weights_RDP, RDP_gamma, imgsiz(1), imgsiz(2), Ndx, Ndy, neighbour_offsets);
    grad = double(grad);
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
function grad = RDP_simulation(x, weights_RDP, gamma, Nx, Ny, Ndx, Ndy, neighbour_offsets)
x = padding_simulation(reshape(x,Nx,Ny),[Ndx Ndy]);
x = x(:);
weights = [weights_RDP(1:ceil(numel(weights_RDP(:))/2) - 1);weights_RDP(ceil(numel(weights_RDP(:))/2) + 1 : end)];
temp1 = x(neighbour_offsets(:,ceil(numel(weights_RDP(:))/2)));
rjk = single(zeros(size(neighbour_offsets, 1), size(neighbour_offsets, 2)));
all_member = [(1:(ceil(numel(weights_RDP(:))/2)-1)) ((ceil(numel(weights_RDP(:))/2)+1) : (size(neighbour_offsets, 2)))];
for i = all_member
    temp2 = x(neighbour_offsets(:, i));
    rjk(:, i) = temp1 ./ temp2;
end
clear neighbour_offsets temp1 temp2
rjk = rjk(:, all_member);
clear all_member
grad = -(((rjk - 1) .* (gamma .* abs(rjk - 1) + rjk + 3)) ./ (rjk + 1 + gamma .* abs(rjk - 1)).^2) * single(weights);
grad(isinf(grad)) = 0;
grad(isnan(grad)) = 0;
end
function weights_RDP = RDP_weights_simulation(weights)
weights_RDP = weights/sum(weights(~isinf(weights)));
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