function y=p_Kondor(l1,l2,l,pos)
% calculates the bispectrum invariants from the Kondor paper

y=0;

for m=-l:l
    for m1=max([-l1 m-l2]):min([l1 m+l2])
        y=y+ClebschGordan(l1,l2,l,m1,m-m1,m)*conj(fhat(l1,m1,pos))*conj(fhat(l2,m-m1,pos))*fhat(l,m,pos);
    end
end

function f=fhat(l,m,pos)
phi = atan2(pos(:,2),pos(:,1));
r = (pos(:,1).^2+pos(:,2).^2+pos(:,3).^2).^0.5;
theta = acos(pos(:,3)./r);

Y_lm=sqrt((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m))*plm(l,m,cos(theta))'.*exp(sqrt(-1)*m*phi);

%f=sum(exp(-r/0.1).*Y_lm);
f=sum(r.*Y_lm);