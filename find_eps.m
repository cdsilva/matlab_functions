function find_eps(W, epsmin, epsmax)

eps = logspace(epsmin, epsmax, 20);
sumA = zeros(size(eps));

for i=1:length(eps)
    sumA(i) = sum(sum(exp(-W/eps(i))));
end

figure;
loglog(eps, sumA)
