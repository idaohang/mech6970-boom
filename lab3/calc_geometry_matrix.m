function G = calc_geometry_matrix(user_pos, svpos)

% CALC_GEOMETRY_MATRIX computes the geometry matrix used for least squares
% estimation of receiver position from pseudorange data
% 
% USAGE:
%   G = calc_geometry_matrix(user_pos, svpos)
% 
% INPUTS:
% 
%   user_pos: ECEF position estimate of the user
%     (1x3 vector)
% 
%   svpos: ECEF position estimates of the satellites
%     (mx3 matrix, where m is the number of sv's)
% 
% OUTPUTS:
% 
%   G: unit vectors from user to sv
%     (mx3 matrix, where m is the number of sv's)
% 
% @author Robert Cofield
% 

sz = size(svpos);
m = sz(1);
if (sz(2)~=3), error('input dimensions'); end
if (any(size(user_pos)~=[1 3])), user_pos = user_pos'; end
  
G = zeros(m,3);
for k = 1:m
  dpos = svpos(k,:) - user_pos;
  G(k,:) = -dpos/norm(dpos,2);
end
