function x = SART_simulation(G, x0, imgsiz, voxel_size, prjsiz, yi, ni, ri, maxit, use_psf, FWHM)
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
wx   = sens;
nn = 2 ^ ceil(log2(prjsiz(2)));
ii = bin2dec(fliplr(dec2bin((0 : (nn - 1)))));
ii = ii(ii < prjsiz(2));
proj_angle_start = 1 + ii;
clear nn ii
y_norm = proj_forw_simulation(G, ones(size(x)), imgsiz, voxel_size, use_psf, FWHM);
for it = 1 : maxit
    disp(['iteration: ' num2str(it)]);
    for i_angle = 1 : prjsiz(2)
        i_proj_angle = (proj_angle_start(i_angle) : prjsiz(2) : prjsiz(2))';
        i_proj_bin = zeros(length(i_proj_angle) * prjsiz(1), 1);
        for i_length_proj_angle = 1 : length(i_proj_angle)
            i_proj_bin(((i_length_proj_angle - 1) * prjsiz(1) + 1) : (i_length_proj_angle * prjsiz(1))) = ...
                ((i_proj_angle(i_length_proj_angle) - 1) * prjsiz(1) + 1) : 1 : (i_proj_angle(i_length_proj_angle) * prjsiz(1));
        end
        yi_proj_angle = yi(i_proj_bin);
        ni_proj_angle = ni(i_proj_bin);
        ri_proj_angle = ri(i_proj_bin);
        y_norm_proj_angle = y_norm(i_proj_bin);
        G_proj_angle = G(i_proj_bin, :);
        yb_proj_angle = ni_proj_angle .* proj_forw_simulation(G_proj_angle, x, imgsiz, voxel_size, use_psf, FWHM) + ri_proj_angle;
        yy_proj_angle_diff = (yi_proj_angle - yb_proj_angle) ./ y_norm_proj_angle;
        xb = proj_back_simulation(G_proj_angle, ni_proj_angle.*yy_proj_angle_diff, imgsiz, voxel_size, use_psf, FWHM);
        x = x + xb ./ wx;
        x(~mask) = 0;
    end
end