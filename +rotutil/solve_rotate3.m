function varargout = solve_rotate3(R)
%% SOLVE_ROTATE3 solves a 3D rotation matrix for its 3 constituent Euler angles
% 
% https://truesculpt.googlecode.com/hg-history/38000e9dfece971460473d5788c235fbbe82f31b/Doc/rotation_matrix_to_euler.pdf
% 
% This assumes that the rotation matrix was formed by rotating about the
% axes in x,y,z order. It presents the two possible solutions for the
% matrix. The first solution represents a y rotation that lies in either
% the first or second quandrant. The second solution represents a y
% rotation that lies in either the third or fourth quadrant.
% 
% USAGE:
%   [x,y,z] = solve_rotate3(R)
%   [xyz1, xyz2] = solve_rotate3(R)
%   [x1, x2, y1, y2, z1, z2] = solve_rotate3(R)
% 
% INPUTS:
%   R: 3x3 x,y,z rotation matrix
% 
% OUTPUTS:
%   x,y,z: radian angles
%
%
if all(R(3,1)~=[-1,1])
  y1 = -asin(R(3,1));
  x1 = atan2( R(3,2)/cos(y1) , R(3,3)/cos(y1) );
  z1 = atan2( R(2,1)/cos(y1) , R(1,1)/cos(y1) );
  
  y2 = pi - y1;
  x2 = atan2( R(3,2)/cos(y2) , R(3,3)/cos(y2) );
  z2 = atan2( R(2,1)/cos(y2) , R(1,1)/cos(y2) );
  x = x2; y = y2; z = z2;
else
  z = 0; % z = anything really.
  if R(3,1) == -1
    y = pi/2;
    x = z + atan2(R(1,2),R(1,3));
  else % R(3,1) == 1
    y = -pi/2;
    x = -z + atan2(-R(1,2),-R(1,3));
  end
  z1=z; z2=z; y1=y; y2=y; x1=x; x2=x;  
end

if nargout == 2
  varargout{1} = [x1 y1 z1]; varargout{2} = [x2 y2 z2];
elseif nargout == 3
  varargout{1} = [x1 x2]; varargout{2} = [y1 y2]; varargout{3} = [z1 z2];
elseif nargout == 6
  varargout{1} = [x1 x2 y1 y2 z1 z2];
else
  error('Incorrect number of output arguments.');
end

end