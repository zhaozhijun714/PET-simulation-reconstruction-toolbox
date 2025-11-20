function tr_offsets = compute_offsets_simulation(image_size, Ndx, Ndy)
image_size = double(image_size);
N = image_size(1) * image_size(2);
s = [image_size(1) + Ndx * 2 image_size(2) + Ndy * 2];
N_pad = min(2, Ndx + Ndy);
[c1{1 : N_pad}] = ndgrid(1 : (Ndx * 2 + 1));
c2(1 : N_pad) = {Ndy + 1};
offsets = sub2ind(s, c1{:}) - sub2ind(s, c2{:});
tr_ind = sub2ind([image_size(1)+Ndx*2 image_size(2)+Ndy*2], ...
    mod((1:N)'-1, image_size(1))+(Ndx+1), mod(floor(((1:double(N))'-1)/image_size(1)), image_size(2))+(Ndy+1));
tr_offsets = uint32(bsxfun(@plus, tr_ind, offsets(:)'));
end

