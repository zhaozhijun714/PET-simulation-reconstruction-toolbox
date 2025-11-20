function x = ISRA_simulation(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM)
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
temp = proj_back_simulation(G, ni.*yi, imgsiz, voxel_size, use_psf, FWHM);
for it = 1 : maxit
    disp(sprintf('iteration %d',it));
    yb = ni .* proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM) + ri;
    xb = proj_back_simulation(G, ni.*yb, imgsiz, voxel_size, use_psf, FWHM);
    x  = x .* (temp ./ xb);
    x(~mask) = 0;
end