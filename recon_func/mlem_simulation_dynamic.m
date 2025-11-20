function x = mlem_simulation_dynamic(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM)
numpix = prod(imgsiz);
numfrm = size(yi,2);
if isempty(x0)
    x0 = ones(numpix,numfrm);
end
for i_frame = 1 : numfrm
    sens(:, i_frame) = proj_back_simulation(G, ni(:, i_frame), imgsiz, voxel_size, use_psf, FWHM);
end
mask = sens > 0;
if isempty(maxit)
    maxit = 10;
end
x    = max(mean(x0(:))*1e-9,x0); x(~mask) = 0;
yeps = mean(yi(:))*1e-9;
wx   = sens;
for it = 1:maxit
    disp(sprintf('iteration %d',it));
    for i_frame = 1 : numfrm
        yb(:, i_frame) = ni(:, i_frame) .* proj_forw_simulation(G, x(:, i_frame), imgsiz, voxel_size, use_psf, FWHM) + ri(:, i_frame);
    end
    yy = yi ./ (yb + yeps);
    yy(yb == 0 & yi == 0) = 1;
    for i_frame = 1 : numfrm
        xb(:, i_frame) = proj_back_simulation(G, ni(:, i_frame) .* yy(:, i_frame), imgsiz, voxel_size, use_psf, FWHM);
    end
    x  = x ./ wx .* xb;
    x(~mask) = 0;
end