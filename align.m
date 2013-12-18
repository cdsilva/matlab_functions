function theta = align(data)

[m, n] = size(data);

theta_ij = zeros(m);

for i=1:m
    for j=1:i-1
        sse = inf;
        for k=1:n
            temp_sse = sum((data(i,:) - circshift(data(j,:),[0 k])).^2);
            if temp_sse < sse
                sse = temp_sse;
                theta_ij(i,j) = k;
                theta_ij(j,i) = -k;
            end
        end
    end
end

theta_ij = 2*pi*theta_ij/n;

% for i=1:m
%     sse(i,i) = 0;
% end

H = exp(-sqrt(-1)*theta_ij);
% H = H .* exp(-sse);

[z, ~] = eigs(H, 1);

z = z ./ abs(z);

theta = log(z) / sqrt(-1) * (n / (2*pi));

theta = round(real(theta));

