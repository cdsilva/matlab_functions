function [V, D, a] = align_dmaps(data, eps)

[m, n] = size(data);

a_ij = zeros(m);

all_rot = zeros(n);
W = zeros(m);

for i=1:m
    all_rot = gallery('circul',data(i,:));
%     for k=1:n
%         all_rot(k,:) = circshift(data(i,:),[0 k]);
%     end
    for j=1:i-1
        [min_dist, ind] = min(sum((all_rot - ones(n,1)*data(j,:)).^2, 2));
        ind = ind - 1;
        a_ij(i, j) = ind;
        a_ij(j, i) = -ind;
        W(i,j) = min_dist;
        W(j,i) = min_dist;
    end
end

T = exp(-2 * sqrt(-1) * pi * a_ij / n);
K = exp(-W/eps);
for i=1:m
    K(i,:) = K(i,:) / sum(K(i,:));
end

S = K.*T;

[V, D] = eigs(S, 10);
[~, I] = sort(abs(diag(D)), 'descend');
V = V(:,I);
D = D(I,I);

a = round(n/(2*pi) * atan2(imag(V(:,1)),real(V(:,1))));

for i=1:m
    V(i,:) = V(i,:) / V(i,1);
end


