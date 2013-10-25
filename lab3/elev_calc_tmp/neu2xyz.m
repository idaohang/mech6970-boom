function [XYZ,CXYZ] = neu2xyz(O,V,SV,COR);

% XYZ2NEU	Convert local topocentric into ECEF
%
%		Input:
%		  O = origin vector in ellipsoidal coordinates (lat lon height)
%		  V = position or velocity vector in NEU frame (m or m/yr)
%		  SV = stdev in NEU frame (m or m/yr)
%		  COR = correlations, NE NU EU
%		  (NOTE: O, V, SV, COR can be n x 3 matrices, n = # of sites)
%
%		Output:
%		  XYZ = output in ECEF Cartesian frame (m)
%		  CXYZ = associated covariance (m), format is:
%		         Cxx Cxy Cxz Cyy Cyz Czz
%		  (NOTE: XYZ and CXYZ will be matrices with n rows)
%
%		Call: [XYZ,CXYZ] = neu2xyz(O,V,SV,COR);

% if O is a single point, make it the same size as V
if (size(O,1) == 1)
  lat = ones(size(V,1),1) .* O(1);
  lon = ones(size(V,1),1) .* O(2);
  h = ones(size(V,1),1) .* O(3);
else
  lat = O(:,1); lon = O(:,2); h = O(:,3);
end

% read rest of input
vn = V(:,1); ve = V(:,2); vu = V(:,3);
svn = SV(:,1); sve = SV(:,2); svu = SV(:,3);
cne = COR(:,1); cnu = COR(:,2); ceu = COR(:,3);

% convert position(s) of origin to ECEF
[XR,YR,ZR] = wgs2xyz(lon,lat,h);

% compute sines and cosines
cp = cos(lon.*pi/180); sp = sin(lon.*pi/180); % longitude = phi
cl = cos(lat.*pi/180); sl = sin(lat.*pi/180); % latitude = lam

% for each site
XYZ = [];
CXYZ = [];
for i=1:size(V,1)
  % build the rotation matrix
  R = [ -sl(i)*cp(i)   -sl(i)*sp(i)    cl(i);
         -sp(i)            cp(i)           0;
         cl(i)*cp(i)    cl(i)*sp(i)    sl(i)];

  % apply the rotation
  XYZi = R' * [vn(i);ve(i);vu(i)];

  % svu cannot be zero or R'*CVi*R may return negative variances
  if (svu(i)==0) svu(i)= mean([svn(i) sve(i)]); end;

  % build covariance for that site
  CVi = [svn(i)^2 cne(i)   cnu(i);
         cne(i)   sve(i)^2 ceu(i);
         cnu(i)   ceu(i)   svu(i)^2];

  % propagate covariance
  CXYZi = R' * CVi * R;

  % increment result matrices
  XYZ = [XYZ;
         XYZi(1) XYZi(2) XYZi(3)];
  CXYZ = [CXYZ;
          CXYZi(1,1) CXYZi(1,2) CXYZi(1,3) CXYZi(2,2) CXYZi(2,3) CXYZi(3,3)];
end

