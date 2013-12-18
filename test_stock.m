
window = 10;

n = length(A)-window-1;

W = zeros(n);
f = zeros(n,1);

for i=1:n
    f(i) = A(i+window+1)/A(i+window);
    for j=1:i-1
        W(i,j) = norm(A(i:i+window)-A(j:j+window));
        W(j,i)=W(i,j);
    end
end

eps = 0.1*window;

[V,D] = dmaps_weight(W,1,eps,30);

plot(V(:,2),V(:,3),'.')

figure(2)
plot(diag(D),'.')

figure(3)
plot(V(:,2),f,'.')

W = exp(-W./eps);
u=1:200;
l=200:n;
Wu = W(u,u);
Du=diag(sum(Wu));
fl = f(l);
fu = inv(Du-Wu)*W(u,l)*fl;

%plot(fu,f(u),'.')