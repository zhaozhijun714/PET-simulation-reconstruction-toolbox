function x = mlem_simulation_osl_mrp(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, beta, medx, medy)
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
for it = 1 : maxit
    disp(sprintf('iteration %d',it));
    grad = MRP_2D_simulation(x, medx, medy, imgsiz(1), imgsiz(2));
    wx_osl = wx + beta * grad;
    wx_osl(wx_osl < 0) = 0;
    yb = ni .* proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM) + ri;
    yy = yi ./ (yb + yeps);
    yy(yb == 0 & yi == 0) = 1;
    xb = proj_back_simulation(G, ni.*yy, imgsiz, voxel_size, use_psf, FWHM);
    x  = x ./ wx_osl .* xb;
    x(isnan(x)) = 0;
    x(isinf(x)) = 0;
    x(~mask) = 0;
end
end
function [grad, grad_temp]= MRP_2D_simulation(im, medx, medy, Nx, Ny)
im = im(:);
grad_temp = medfilt2(reshape(im, Nx, Ny), [medx medy], 'symmetric');
grad_temp = grad_temp(:);
grad = (im - grad_temp) ./ grad_temp;
grad = grad(:);
grad(isnan(grad)) = 0;
grad(isinf(grad)) = 0;
end