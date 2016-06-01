function vol = getVolCylRecip(phi,B,S,CR,lambda,angleRD)
if angleRD
    phi = phi/180*pi;
end;
a = S / 2;
l = a / lambda;
aPist = B^2*pi/4;
vDisp = aPist*S;
vComp = vDisp/(CR - 1);
vol = vComp + aPist*(a + l - (a*cos(phi) + sqrt(l^2 - (a*sin(phi)).^2)));
