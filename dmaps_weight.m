function [V,D]=dmaps_weight(W,alpha,eps,n)
%%computes dmap embedding of matrix with specified alpha
%%W_ij=d^2(x_i,x_j) and kernel width
%%eps; returns top n eigenvedctors and eigenvalues sorted according to eigenvalue

W=exp(-W./eps);

d=sum(W);

W=diag(d.^-alpha)*W*diag(d.^-alpha);

d=sum(W);

S=diag(d.^-0.5)*W*diag(d.^-0.5);

if n==max(size(S))
    [V,D]=eig(S);
else
    [V,D]=eigs(S,n);
end

V=diag(d.^-0.5)*V;

[sortd I]=sort(abs(diag(D)),'descend');

D=D(I,I);
V=V(:,I);

for i=1:n
    V(:,i)=V(:,i)./norm(V(:,i));
end



