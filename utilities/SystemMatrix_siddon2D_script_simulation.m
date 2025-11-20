function A_system_matrix = SystemMatrix_siddon2D_script_simulation(scanner, reconstruction)
scanner.Radius = scanner.Radius + reconstruction.DOI;
GridSize = double(reconstruction.image_dim);
voxel_number = prod(GridSize);
VoxelSize = double(reconstruction.voxel_size);
Views = floor(scanner.crystal_per_ring / 2);
S_WIDTH = scanner.crystal_per_ring - scanner.NrModulesTrans * scanner.NrCrystalsTrans; % Number of radial sinogram bins
Geometry = GeometryGenerator_simulation(scanner);
[sinogram_ID1_2D, sinogram_ID2_2D] = calculate_sino_coor_from_id_simulation(scanner, S_WIDTH, Views);
LOR_number = size(sinogram_ID1_2D, 1);
sinoid1_coor = zeros(size(sinogram_ID1_2D, 1), 3, 'double');
sinoid1_coor(:, 1:2) = [Geometry(sinogram_ID1_2D, 1), Geometry(sinogram_ID1_2D, 2)];
sinoid2_coor = zeros(size(sinogram_ID2_2D, 1), 3, 'double');
sinoid2_coor(:, 1:2) = [Geometry(sinogram_ID2_2D, 1), Geometry(sinogram_ID2_2D, 2)];
sinoid_all = (1 : LOR_number)';
voxel_index = cell((LOR_number), 1);
intersection_length = cell((LOR_number), 1);
lor_record = zeros((LOR_number),2,'int32');
for i_lor = 1 : LOR_number
    [voxel_index_temp, intersection_length_temp] = siddon2D_simulation(sinoid1_coor(i_lor, :), sinoid2_coor(i_lor, :), GridSize, VoxelSize);
    voxel_index{i_lor} = int32(voxel_index_temp);
    intersection_length{i_lor} = double(intersection_length_temp);
    if isempty(voxel_index_temp)
        lor_record(i_lor, 1) = 0;
        lor_record(i_lor, 2) = 0;
    else
        lor_record(i_lor, 1) = sinoid_all(i_lor);
        lor_record(i_lor, 2) = length(voxel_index_temp);
    end
end
intersection_length = cell2mat(intersection_length);
voxel_index = cell2mat(voxel_index) - 1;
lor_record = lor_record(:,2);
lor_record = repelem(uint32(1 : length(lor_record)), uint32(lor_record))';
voxel_index = uint32(voxel_index) + 1;
A_system_matrix = sparse(lor_record, voxel_index, double(intersection_length), LOR_number, voxel_number);
end

function [sinogram_ID1_2D, sinogram_ID2_2D] = calculate_sino_coor_from_id_simulation(scanner, Ndist, Nang)
N_DET = scanner.crystal_per_ring;
ID = formDetectorIndices_simulation(N_DET);
ID1 = double(ID(:, 1)) - 1;
ID2 = double(ID(:, 2)) - 1;
clear  ID
sinogram_ID1_2D = zeros(Ndist, Nang, 'single');
sinogram_ID2_2D = zeros(Ndist, Nang, 'single');
for i_ID = 1 : length(ID1)
    crystal1 = mod(ID1(i_ID), N_DET);
    crystal2 = mod(ID2(i_ID), N_DET);
    phi = floor((mod((crystal1 + crystal2 + floor(N_DET / 2)), N_DET)) / 2);
    if (((crystal1 + crystal2) < (floor(3 * N_DET / 2 ))) && ((crystal1 + crystal2) >= floor(N_DET / 2 )))
        u = abs(crystal1 - crystal2) -  floor(N_DET / 2) + floor(Ndist / 2);
    else
        u = -abs(crystal1 - crystal2) +  floor(N_DET / 2) + floor(Ndist / 2);
    end
    if (u >= Ndist) || (u < 0)
        continue
    end
    sinogram_ID1_2D(u + 1, phi + 1) = ID1(i_ID);
    sinogram_ID2_2D(u + 1, phi + 1) = ID2(i_ID);
end
sinogram_ID1_2D = sinogram_ID1_2D(:) + 1;
sinogram_ID2_2D = sinogram_ID2_2D(:) + 1;
end
function L = formDetectorIndices_simulation(det_per_ring)
L = zeros(sum(1 : det_per_ring) - det_per_ring, 2, 'uint32');
idx = 1;
for kk = int32(1) : det_per_ring
    current_indices = [repelem(kk, det_per_ring - kk)', (kk + 1 : det_per_ring)'];
    L(idx:(idx + size(current_indices, 1) - 1), :) = current_indices;
    idx = idx + size(current_indices, 1);
end
L(L(:, 1) == 0, :) = [];
end
function Geometry = GeometryGenerator_simulation(scanner)
NrSectorsTrans = scanner.NrSectorsTrans;
AngleStart = scanner.AngleStart;
Radius = scanner.Radius;
NrModulesTrans = scanner.NrModulesTrans;
ModulesGapTrans = scanner.ModulesGapTrans;
NrCrystalsTrans = scanner.NrCrystalsTrans;
DetectorSize = scanner.CrystalSize;
DetectorGap = scanner.CrystalGap;
NrCrystalsPerSectorTrans = NrModulesTrans * NrCrystalsTrans;
NrCrystalsPerRing = NrCrystalsPerSectorTrans * NrSectorsTrans;
AngleStep = 360.0 / NrSectorsTrans;
Geometry = zeros(NrCrystalsPerRing,2, 'single');
for Sector = 1:NrSectorsTrans
    Angle = 90 + AngleStart + (Sector-1) * AngleStep;
    AngleSector = Angle - 90;
    Xmid = Radius * cosd(Angle);
    Ymid = Radius * sind(Angle);
    if mod(NrModulesTrans,2) == 0
        Shift = .5*(ModulesGapTrans + DetectorSize(1));
    else
        Shift = .5 * (DetectorSize(1) + DetectorGap(1));
    end
    Xdetector = Xmid + Shift * cosd(AngleSector);
    Ydetector = Ymid + Shift * sind(AngleSector);
    Detector = floor((Sector - .5) * NrCrystalsPerSectorTrans + 1);
    while Detector < (Sector * NrCrystalsPerSectorTrans + 1)
        Geometry(Detector,1) = Xdetector;
        Geometry(Detector,2) = Ydetector;
        if round(Detector / NrCrystalsTrans) == (Detector / NrCrystalsTrans)
            Xdetector = Xdetector + (DetectorSize(1) + ModulesGapTrans) * cosd(AngleSector);
            Ydetector = Ydetector + (DetectorSize(1) + ModulesGapTrans) * sind(AngleSector);
        else
            Xdetector = Xdetector + (DetectorSize(1) + DetectorGap(1)) * cosd(AngleSector);
            Ydetector = Ydetector + (DetectorSize(1) + DetectorGap(1)) * sind(AngleSector);
        end
        Detector = Detector + 1;
    end
    Xdetector = Xmid - Shift * cosd(AngleSector);
    Ydetector = Ymid - Shift * sind(AngleSector);
    Detector = floor((Sector - .5) * NrCrystalsPerSectorTrans);
    while Detector > (Sector - 1) * NrCrystalsPerSectorTrans
        Geometry(Detector,1) = Xdetector;
        Geometry(Detector,2) = Ydetector;
        if round(((Detector - 1) / NrCrystalsTrans)) == ((Detector - 1) / NrCrystalsTrans)
            Xdetector = Xdetector - (DetectorSize(1) + ModulesGapTrans) * cosd(AngleSector);
            Ydetector = Ydetector - (DetectorSize(1) + ModulesGapTrans) * sind(AngleSector);
        else
            Xdetector = Xdetector - (DetectorSize(1) + DetectorGap(1)) * cosd(AngleSector);
            Ydetector = Ydetector - (DetectorSize(1) + DetectorGap(1)) * sind(AngleSector);
        end
        Detector = Detector - 1;
    end
    Geometry(Detector + 1 : Detector + NrCrystalsPerSectorTrans, :) = flip(squeeze(Geometry( Detector + 1 : Detector + NrCrystalsPerSectorTrans, :)), 1);
end
if strcmp(scanner.Generate_crystal_direction, 'CW')
    Geometry(: , :) = flip(squeeze(Geometry(:, :)), 1);
end
end
function [num_voxel, len] = siddon2D_simulation(xstart, xend, image_dim, voxel_size)
nx = image_dim(1);
ny = image_dim(2);
dx = voxel_size(1);
dy = voxel_size(2);
Nx = nx + 1;
Ny = ny + 1;
Xplane1 = -nx * dx / 2;
Yplane1 = -ny * dy / 2;
Xplanen = nx * dx / 2;
Yplanen = ny * dy / 2;
alphax1 = (Xplane1 - xstart(1)) / (xend(1) - xstart(1));
alphaxn = (Xplanen - xstart(1)) / (xend(1) - xstart(1));
alphay1 = (Yplane1 - xstart(2)) / (xend(2) - xstart(2));
alphayn = (Yplanen - xstart(2)) / (xend(2) - xstart(2));
alphamin = max(max(0, min(alphax1, alphaxn)), max(min(alphay1, alphayn)));
alphamax = min(min(1, max(alphax1, alphaxn)), min(max(alphay1, alphayn)));
if alphamax >= alphamin
    if xend(1) >= xstart(1)
        lmin = ceil(Nx - (Xplanen - alphamin * (xend(1) - xstart(1)) - xstart(1)) / dx);
        lmax = floor(1 + (xstart(1) + alphamax * (xend(1) - xstart(1)) - Xplane1) / dx);
        alphax = zeros(1, lmax - lmin + 1);
        for l = 1 : lmax - lmin + 1
            alphax(l) = (Xplane1 + (l + lmin - 2) * dx - xstart(1)) / (xend(1) - xstart(1));
        end
    else
        lmin = ceil(Nx - (Xplanen - alphamax * (xend(1) - xstart(1)) - xstart(1)) / dx);
        lmax = floor(1 + (xstart(1) + alphamin * (xend(1) - xstart(1)) - Xplane1) / dx);
        alphax = zeros(1, lmax - lmin + 1);
        for l = lmax - lmin + 1 : -1 : 1
            alphax(lmax - lmin + 2 - l) = (Xplane1 + (l + lmin - 2) * dx - xstart(1)) / (xend(1) - xstart(1));
        end
    end
    if xend(2) >= xstart(2)
        mmin = ceil(Ny - (Yplanen - alphamin * (xend(2) - xstart(2)) - xstart(2)) / dy);
        mmax = floor(1 + (xstart(2) + alphamax * (xend(2) - xstart(2)) - Yplane1) / dy);
        alphay = zeros(1, mmax - mmin + 1);
        for m = 1 : mmax - mmin + 1
            alphay(m) = (Yplane1 + (m + mmin - 2) * dy - xstart(2)) / (xend(2) - xstart(2));
        end
    else
        mmin = ceil(Ny - (Yplanen - alphamax * (xend(2) - xstart(2)) - xstart(2)) / dy);
        mmax = floor(1 + (xstart(2) + alphamin * (xend(2) - xstart(2)) - Yplane1) / dy);
        alphay = zeros(1, mmax - mmin + 1);
        for m = mmax - mmin + 1 : -1 : 1
            alphay(mmax - mmin + 2 - m) = (Yplane1 + (m + mmin - 2) * dy - xstart(2)) / (xend(2) - xstart(2));
        end
    end
    alpha_temp = [alphax, alphay];
    alpha = sort(alpha_temp);
    d = sqrt((xstart(1) - xend(1)) ^ 2 + (xstart(2) - xend(2)) ^ 2);
    l1 = alpha(1 : size(alpha, 2) - 1);
    l2 = alpha(2 : size(alpha, 2));
    len = d * (l2 - l1);
    voxel_i = floor(1 + (xstart(1) + (l1 + l2) * (xend(1) - xstart(1)) / 2 - Xplane1) / dx) - 1;
    voxel_j = floor(1 + (xstart(2) + (l1 + l2) * (xend(2) - xstart(2)) / 2 - Yplane1) / dy) - 1;
    num_voxel = voxel_j * nx + voxel_i + 1;
    len = len';
    num_voxel = num_voxel';
else
    num_voxel = [];
    len = [];
end
end