function y = proj_forw_simulation(G, x, imgsiz, voxel_size, use_psf, FWHM)
if use_psf
    x_psf = PSFmodelFilter_simulation(x, imgsiz, voxel_size, FWHM);
else
    x_psf = x;
end
y = G * x_psf;
end