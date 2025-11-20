function weights = compute_weights_simulation(voxel_size, Ndx, Ndy)
weights = zeros(((Ndx * 2 + 1) * (Ndy * 2 + 1)), 1);
edist = zeros((Ndx * 2 + 1), 1);
cc = zeros((Ndx * 2+ 1) * (Ndy * 2 + 1),1);
ll = 0;
for kk = Ndy : -1 : -Ndy
    ll = ll + 1;
    apu = [((Ndx : -1 : -Ndx) * voxel_size(1))', (repelem(kk, Ndy * 2 + 1) * voxel_size(2))'];
    for ii = 1 : length(apu)
        edist(ii) = sqrt(apu(ii, :) * apu(ii, :)');
    end
    cc((Ndy * 2 + 1) * (ll - 1) + 1 : (Ndy * 2 + 1) * ll) = edist;
end
weights(1: (Ndx * 2 + 1) * (Ndy * 2 + 1)) = cc;
weights = 1 ./ weights;
end

