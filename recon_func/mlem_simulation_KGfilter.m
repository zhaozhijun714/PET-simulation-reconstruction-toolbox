function x = mlem_simulation_KGfilter(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, kernel_filter, use_psf, FWHM)
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end
if isempty(kernel_filter)
    kernel_filter = speye(numpix);
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
    yb = ni .* proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM) + ri;
    yy = yi ./ (yb + yeps);
    yy(yb == 0 & yi == 0) = 1;
    xb = proj_back_simulation(G, ni.*yy, imgsiz, voxel_size, use_psf, FWHM);
    x  = x ./ wx .* xb;
    x(~mask) = 0;
    x = kernel_filter * x;
end