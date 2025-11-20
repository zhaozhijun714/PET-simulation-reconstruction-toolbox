function xb = proj_back_simulation(G, yi, imgsiz, voxel_size, use_psf, FWHM)
xb = G' * yi;
if use_psf
    xb = PSFmodelFilter_simulation(xb, imgsiz, voxel_size, FWHM);
end
end