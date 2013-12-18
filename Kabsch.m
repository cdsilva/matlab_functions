function Pnew = Kabsch(P, Q)

A = P'*Q;

[V, S, W] = svd(A);

d = sign(det(W*V'));

I = eye(size(P,2));
I(end,end) = d;

U = W*I*V';

Pnew = P*U';
