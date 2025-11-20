function x = osem_simulation(G, x0, imgsiz, voxel_size, prjsiz, yi, ni, ri, maxit, subset_number, use_psf, FWHM)
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
nn = 2 ^ ceil(log2(subset_number));
ii = bin2dec(fliplr(dec2bin((0 : (nn - 1)))));
ii = ii(ii < subset_number);
subset_start = 1 + ii;
clear nn ii
for it = 1 : maxit
    for i_sub = 1 : subset_number
        iter = it * (subset_number - 1) + i_sub;
        disp(['iteration: ' num2str(it) ', subset: ' num2str(i_sub)]);
        i_subset_angle = (subset_start(i_sub) : subset_number : prjsiz(2))';
        i_subset_bin = zeros(length(i_subset_angle) * prjsiz(1), 1);
        for i_length_subset_angle = 1 : length(i_subset_angle)
            i_subset_bin(((i_length_subset_angle - 1) * prjsiz(1) + 1) : (i_length_subset_angle * prjsiz(1))) = ...
                ((i_subset_angle(i_length_subset_angle) - 1) * prjsiz(1) + 1) : 1 : (i_subset_angle(i_length_subset_angle) * prjsiz(1));
        end
        yi_subset = yi(i_subset_bin);
        ni_subset = ni(i_subset_bin);
        ri_subset = ri(i_subset_bin);
        G_subset = G(i_subset_bin, :);
        wx_subset = full(proj_back_simulation(G_subset, ni_subset, imgsiz, voxel_size, use_psf, FWHM));
        yb_subset = ni_subset .* proj_forw_simulation(G_subset, x, imgsiz, voxel_size, use_psf, FWHM) + ri_subset;
        yy_subset = yi_subset ./ (yb_subset + yeps);
        yy_subset(yb_subset == 0 & yi_subset == 0) = 1;
        xb = proj_back_simulation(G_subset, ni_subset.*yy_subset, imgsiz, voxel_size, use_psf, FWHM);
        x  = x ./ wx_subset .* xb;
        x(~mask) = 0;
    end
end