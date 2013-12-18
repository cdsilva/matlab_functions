function [V,D]=dmaps(W,eps,n,Wtol)
% computes dmap embedding of matrix with W_ij=d^2(x_i,x_j) 
% kernel width eps 
% returns top n eigenvedctors and eigenvalues sorted by eigenvalue
% Wtol is an optional parameter; if Wtol is set, then entries in the kernel
% matrix that are < Wtol will be set to 0 and the matrix will be treated as
% a sparse matrix; this allows for larger matrices
    
W=exp(-W./eps);

if nargin > 3
    ind = W < Wtol;
    W(ind) = 0;
end

d=sum(W);

S=diag(d.^-0.5)*W*diag(d.^-0.5);

if nargin > 3
    S = sparse(S);
end

if n==max(size(S))
    [V,D]=eig(S);
else
    %opts.maxit = 3000;
    %opts.tol = 1e-8;
    OPTS.issym = true;
    [V,D]=eigs(S,n,'LM',OPTS);
end

V=diag(d.^-0.5)*V;

[~,I]=sort(abs(diag(D)),'descend');

D=D(I,I);
V=V(:,I);

for i=1:n
    V(:,i)=V(:,i)./norm(V(:,i));
end

