function [V,D]=PCA(data,n_evecs)
%%computes PCs of data

[m, n] = size(data);

for i=1:n
    data(:,i) = data(:,i) - mean(data(:,i));
end

CV = data' * data;

if n_evecs == n
    [V, D] = eig(CV);
else
    [V, D] = eigs(CV, n_evecs);
end

[~,I]=sort(abs(diag(D)),'descend');

D=D(I,I);
V=V(:,I);

for i=1:n_evecs
    V(:,i)=V(:,i)./norm(V(:,i));
end

