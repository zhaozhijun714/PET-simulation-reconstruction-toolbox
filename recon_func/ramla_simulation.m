function x = ramla_simulation(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, relaxation_parameter)
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end
if isempty(maxit)
    maxit = 10;
end
sens = proj_back_simulation(G, ni, imgsiz, voxel_size, use_psf, FWHM);
mask = sens > 0;
x    = max(mean(x0(:))*1e-9,x0(:)); x(~mask) = 0;
yeps = mean(yi(:))*1e-9;
wx   = sens;
lam = zeros(maxit, 1);
lam(1) = relaxation_parameter;
for i = 1 : maxit
    lam(i + 1) = lam(1) / ((i - 1)/20 + 1);
end
lam = lam .* (1 ./ max(wx, [], 'all'));
for it = 1:maxit
    disp(sprintf('iteration %d',it));
    yb = ni .* proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM) + ri;
    yy = yi ./ (yb + yeps) - 1;
    yy(yb == 0 & yi == 0) = 1;
    xb = proj_back_simulation(G, ni.*yy, imgsiz, voxel_size, use_psf, FWHM);
    x = x + lam(it) .* x .* xb;
    x(~mask) = 0;
    x(x < 0) = 0;
end