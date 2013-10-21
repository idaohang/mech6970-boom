function varargout = calc_azel(userpos, svpos)

% CALC_AZEL finds the azimuth and elevation of each satellite given a user
% position.
% USAGE:
%   ae = CALC_AZEL(userpos, svpos)
%   [az,el] = CALC_AZEL(userpos, svpos)
% INPTUTS:
%   userpos: 1x3 vector of ECEF (x,y,z) user position
%   svpos: Nx3 matrix of ECEF(x,y,z) sv positions
%     rows correspond the individual SV's
% OUTPUTS:
%   - All output angles are in radians
%   ae: [azi1  azi2 ... aziN]
%       [ele1  ele2 ... eleN]
%   az: 1st row of ae
%   el: 2nd row of ae
% DEPENDENCIES:
% 

nsv = size(svpos); nsv = nsv(1);

rays_ecef = svpos - repmat(userpos, nsv,1);
[reflat, reflon, ~] = wgsxyz2lla(userpos, 100);
az = zeros(1,nsv);
el = zeros(1,nsv);
for k = 1:nsv
  enu = rotxyz2enu(rays_ecef(k,:)', reflat,reflon);
  az(k) = atan2(enu(1),enu(2));
%   gnd_proj = norm(enu(1:2),2);
  svdist = norm(enu);
  el(k) = asin(enu(3)/svdist);
end

if nargout == 1
  varargout{1} = [az;el];
elseif nargout == 2
  varargout{1} = az;
  varargout{2} = el;
else
  error('Output numbers');
end

end