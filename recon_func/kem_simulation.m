function x = kem_simulation(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, K, use_psf, FWHM)
numpix = prod(imgsiz);
if isempty(x0)
    x0 = ones(numpix,1);
end
if isempty(K)
    K = speye(numpix);
end
if isempty(maxit)
    maxit = 10;
end
sens = proj_back_simulation(G, ni, imgsiz, voxel_size, use_psf, FWHM);
sens = ker_back_simulation(K, sens);
mask = sens > 0;
x    = max(mean(x0(:)) * 1e-9, x0(:));
x(~mask) = 0;
yeps = mean(yi(:)) * 1e-9;
wx   = sens;
for it = 1 : maxit
    disp(sprintf('iteration %d',it));
    z  = ker_forw_simulation(K, x);
    yb = ni .* proj_forw_simulation(G, z, imgsiz, voxel_size, use_psf, FWHM) + ri;
    yy = yi ./ (yb + yeps);
    yy(yb == 0 & yi == 0) = 1;
    zb = proj_back_simulation(G, ni.*yy, imgsiz, voxel_size, use_psf, FWHM);
    xb = ker_back_simulation(K, zb);
    x  = x ./ wx .* xb;
    x(~mask) = 0;
end
end
function y = ker_forw_simulation(K, x)
y = K * x;
end
function x = ker_back_simulation(K, y)
x = K' * y;
end


    
