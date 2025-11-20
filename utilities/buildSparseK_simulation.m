function K = buildSparseK_simulation(N, W, J, numvox, flag)
% gbwang@ucdavis.edu (05-20-2013)
if nargin<3 | isempty(J)
    J = [1:size(N,1)]';
end
if nargin<4 | isempty(numvox)
    numvox = size(N,1);
end
if nargin<5
    flag = 1;
end
K = sparse(J, J(N(:,1)), W(:,1), numvox, numvox);
for n = 2:size(N,2)
    Kn = sparse(J, J(N(:,n)), W(:,n), numvox, numvox);
    K  = K + Kn;
end
sumK = full(sum(K,2));
I = find(sumK==0);
for i = 1:length(I)
    K(I(i),I(i)) = 1;
end
if flag
    K = ( K' + K )/2;
    K = make_sym_simulation(K);
end
end
function newK = make_sym_simulation(K, maxit)
if nargin<2
    maxit = 100;
end
N = size(K,1);
r = ones(N,1);
for it = 1:maxit
    c = 1./ (K'*r);
    r = 1./ (K *c);
end
C = spdiags(c, 0, N, N);
R = spdiags(r, 0, N, N);
newK = R * K * C;
newK = ( newK' + newK )/2;
end

