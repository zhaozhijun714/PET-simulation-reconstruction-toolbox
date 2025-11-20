
clear; clc;
% image size
imgsiz = [111 111]; 
% Simulate dynamic PET sinograms
% phantom
load("zubal_dynamic_model.mat")
disp('--- Generate system matrix ...')
%%% Information on PET-scanner
scanner.ScannerName = '8panel';
scanner.NrSectorsAxial = 4;          % number of sectors axial - contains modules 
scanner.NrSectorsTrans = 8;         % number of sectors transaxial - contains modules side by side in the same slope angle
scanner.SectorGapAxial = 1;        % sector gap axial
scanner.AngleStart = 0;              % angle of the first sector
scanner.Radius = 143.16;                 % radius PET-scanner
scanner.NrModulesAxial = 3;          % number of modules axial - contains detectors
scanner.NrModulesTrans = 6;          % number of modules transaxial - contains detectors
scanner.ModulesGapAxial = 1;       % module gap axial
scanner.ModulesGapTrans = 1;       % module gap transaxial
scanner.NrCrystalsAxial = 10;         % number of detectors per module axial
scanner.NrCrystalsTrans = 8;         % number of detectors per module transaxial
scanner.CrystalSize = [2.05, 2.05];   % size of detector transaxial and axial
scanner.CrystalGap = [0.2, 0.2];  % detector gap transaxial and axial
scanner.crystal_per_ring            =   scanner.NrCrystalsTrans * scanner.NrModulesTrans * scanner.NrSectorsTrans ;
scanner.crystal_ring                =   scanner.NrCrystalsAxial * scanner.NrModulesAxial * scanner.NrSectorsAxial ;
scanner.crystal_all                 =   scanner.crystal_per_ring * scanner.crystal_ring;
scanner.Generate_crystal_direction = 'CW';
reconstruction.DOI = 0; 
reconstruction.image_dim = int32([imgsiz, 1]); 
reconstruction.voxel_size = [1, 1, 1];
tic
G = SystemMatrix_siddon2D_script_simulation(scanner, reconstruction);
toc
prjsiz = [scanner.crystal_per_ring - scanner.NrModulesTrans * scanner.NrCrystalsTrans, scanner.crystal_per_ring / 2];
% frame durations
dt = (scant(:,2)-scant(:,1))/60;
numfrm = length(dt);
% dynamic images of activity
X0 = zeros(prod(size(model)),numfrm);
ID = [0 32 64 128 196 255];
for m = 1:numfrm
    for n = 1:length(ID)
        X0(model==ID(n),m) = TAC(m,n)*dt(m);
    end
end
% noise-free geometric projection
for m = 1 : size(X0, 2)
    proj(:,m) = G * X0(:,m);
end
% attenuation factor
ai = exp(-repmat(((G * u(:)) ./ 10), [1 numfrm]));  
% background (randoms and scatters)
ri = repmat(mean(ai .* proj, 1) * 0.2, [size(proj, 1) 1]);  % 20% uniform background
% total noiseless projection
y0 = ai.*proj + ri; 
% count level
count = 8e7; % a total of 8 million events
% normalized sinograms
cs = count / sum(y0(:));
y0 = cs * y0;
ri = cs * ri;
yi = poissrnd(y0); % noisy projection
ni = ai*cs; % multiplicative factor
% initial estimate
xinit = [];
maxit = 60;
use_psf = 1;
FWHM = 0.75;
voxel_size = double(reconstruction.voxel_size(1:2));
%% Kernelized EM 
disp('--- Building image prior and kernel matrix ...')
% build image prior using 3 composite images
M = {[1:16], [17:20], [21:24]};
for i = 1:length(M)
    y_i = sum(yi(:,M{i}),2);
    n_i = sum(ni(:,M{i}),2);
    r_i = sum(ri(:,M{i}),2);
    x_i = mlem_simulation(G, xinit, imgsiz, voxel_size, y_i, n_i, r_i, 50, 0, FWHM); 
    x_i = gaussfilt_simulation(x_i,imgsiz,3);
    U(:,i) = x_i(:) / sum(dt(M{i}));
end
U = U * diag(1./std(U,1)); % normalization
% build the kernel matrix K using k nearest neighbors
tic;
sigma = 1;
[N, W] = buildKernel_simulation(imgsiz, 'knn', 48, U, 'radial', sigma);
Ks = buildSparseK_simulation(N, W);
toc;    
% shift-invariant temporal kernel
wsize = 3;
[N1, W1] = buildKernelT_simulation(numfrm, 'cube', wsize, [], 'radial');
Kt = buildSparseK_simulation(N1, W1);
% % data-driven temporal kernel
% wsize = 3;
% U = ((yi-ri)./ni);
% for m = 1:numfrm
%     u = gaussfilt(U(:,m),prjsiz,7);
%     U(:,m) = u(:);
% end
% x = mean(U,2); mask = x>mean(x);
% [N1, W1] = buildKernelT_New(numfrm, 'cube', wsize, U(mask,:)', 'radial', 1.0, 1, 1);
% Kt = buildSparseK(N1, W1);
disp('--- STKEM reconstruction ...')
a = stkem_simulation_dynamic(G, xinit, imgsiz, voxel_size, yi, ni, ri, maxit, use_psf, FWHM, Ks, Kt);
x = Ks * a * Kt';
for i_frame = 1 : numfrm
    figure
    imagesc(reshape(x(:, i_frame), imgsiz))
    axis image; axis off
    title(['STKEM Frame : ' num2str(i_frame)])
end