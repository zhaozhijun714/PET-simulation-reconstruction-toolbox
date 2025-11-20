
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
maxit = 1000;
recon_frame = [2, 12, 24];
use_psf = 1;
FWHM = 0.75;
voxel_size = double(reconstruction.voxel_size(1:2));
disp('--- SIRT reconstruction ...')
for i_reocn_frame = 1 : length(recon_frame)
    recon_frame_temp = recon_frame(i_reocn_frame);
    x(:, i_reocn_frame) = SIRT_simulation(G, xinit, imgsiz, voxel_size, yi(:, recon_frame_temp), ni(:, recon_frame_temp), ri(:, recon_frame_temp), maxit, use_psf, FWHM);
end
for i_frame = 1 : length(recon_frame)
    figure
    imagesc(reshape(x(:, i_frame), imgsiz))
    axis image; axis off
    title(['SIRT Frame : ' num2str(recon_frame(i_frame))])
end