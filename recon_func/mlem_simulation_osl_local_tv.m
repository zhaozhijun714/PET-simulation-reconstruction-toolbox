function x = mlem_simulation_osl_local_tv(G, x0, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, beta, niter_tv, gradient_descent_step)
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
    grad = local_TV2D_grad_simulation(reshape(x, imgsiz), gradient_descent_step, niter_tv);
    grad = grad(:);
    grad(isinf(grad)) = 0;
    grad(isnan(grad)) = 0;
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
function f = local_TV2D_grad_simulation(f, dtvg, niter)
for ii = 1:niter
    df = gradientTVnormForward(f);
    df = df ./ sqrt(mean(df(:).^2));
    f = f - dtvg .* df;
end
end
function tvg = gradientTVnormForward(f)
Gx = diff(f, 1, 1);
Gy = diff(f, 1, 2);
Gx = cat(1, Gx, zeros(1, size(Gx, 2), class(f)));
Gy = cat(2, Gy, zeros(size(Gy, 1), 1, class(f)));
nrm = sqrt(Gx.^2 + Gy.^2) + 1e-7;
tvg = Gx([1, 1:end-1], :) - Gx + Gy(:, [1, 1:end-1]) - Gy;
tvg = tvg ./ nrm;
end

