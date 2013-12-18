function R_opt = ang_synch(R, dim)

n = size(R,1)/dim;

[V, D] = eigs(R, dim);

R_opt = zeros(n*dim, dim);

for i=1:n
    [u, s, v] = svd(V(dim*(i-1)+1:dim*i,1:dim));
    R_opt(dim*(i-1)+1:dim*i,:) = u*v';
end
