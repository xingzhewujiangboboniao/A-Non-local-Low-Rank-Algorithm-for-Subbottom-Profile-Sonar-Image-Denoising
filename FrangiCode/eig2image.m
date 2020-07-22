function [Lambda1,Lambda2,Ix,Iy]=eig2image(Dxx,Dxy,Dyy)

tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
v2x = 2*Dxy;
v2y = Dyy - Dxx + tmp;

mag = sqrt(v2x.^2 + v2y.^2); 
i = (mag ~= 0);
v2x(i) = v2x(i)./mag(i);
v2y(i) = v2y(i)./mag(i);

v1x = -v2y; 
v1y = v2x;

mu1 = 0.5*(Dxx + Dyy + tmp);
mu2 = 0.5*(Dxx + Dyy - tmp);

check=abs(mu1)>abs(mu2);

Lambda1=mu1; 
Lambda1(check)=mu2(check);
Lambda2=mu2; 
Lambda2(check)=mu1(check);

Ix=v1x; 
Ix(check)=v2x(check);
Iy=v1y; 
Iy(check)=v2y(check);
end



