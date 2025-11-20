function [image_add_mask, mask] = image_mask_simulation(image_raw, img_dim, voxel_size, fov_transaxial)
Rad = fov_transaxial * 0.5;
img_dim = single(img_dim);
image_raw = reshape(image_raw, img_dim);
[XX, YY] = meshgrid((1 : img_dim(2)) - (img_dim(2) + 1) / 2, (1 : img_dim(1)) - (img_dim(1) + 1) / 2);
XX = XX * voxel_size(1);
YY = YY * voxel_size(2);
XX = XX .^ 2;
YY = YY .^ 2;
R2 = Rad .^ 2;
mask_xy = (XX + YY) - R2;
mask_xy(mask_xy > 0) = 0;
mask_xy(mask_xy < 0) = 1;
mask = mask_xy > 0;
mask = double(mask);
image_add_mask = image_raw .* mask;
image_add_mask = image_add_mask(:);
end