function Ro = interp_3Drotmat(xi, Ri, xo, method, varargin)
% Interpolates a 3x3xn array of 3x3 rotation matrices from xi to xo
% 
% USAGE:
%   Ro = interp_3Drotmat(xi, Ri, xo, method)
%   Ro = interp_3Drotmat(xi, Ri, xo, method, 2)
%     - this uses the second solution to all rotation matrices
% 
% INPUTS:
%   xi:
%   Ri:
%   xo:
%   method:
% 
% OUTPUTS:
%   Ro:
%

if nargin==5
  rot_soln = varargin{1};
  spec_rot_soln = true;
else
  spec_rot_soln = false;
end

leni = length(xi);
leno = length(xo);
Ro = zeros(3,3,leno);
ei = zeros(3,leni); % euler angles for input array
eo = zeros(3,leno); % euler angles for output array

% determine input euler angles
for n = 1:leni
  [x,y,z] = rotutil.solve_rotate3(Ri(:,:,n)); % will return 2 solution sets
  % use the euler angle solution which has the smallest change between
  % epochs
  if n==1 && ~spec_rot_soln
    x = x(1); y = y(1); z = z(1); % this could probably be done better
  elseif spec_rot_soln
    x = x(rot_soln); y = y(rot_soln); z = z(rot_soln);
  elseif ~spec_rot_soln
    [m, rot_soln] = min(abs( (z+2*pi)-(ei(3,n-1)+2*pi) ));
    x=x(rot_soln); y=y(rot_soln); z=z(rot_soln);
  else
    error('?');
  end
  ei(:,n) = [x;y;z];
end

% find output euler angles
eo(1,:) = interp1(xi, ei(1,:), xo, method, 'extrap');
eo(2,:) = interp1(xi, ei(2,:), xo, method, 'extrap');
eo(3,:) = interp1(xi, ei(3,:), xo, method, 'extrap');

% create output rotation matrices
for n = 1:leno
  Ro(:,:,n) = rotutil.rotate3(eo(:,n),'rad');  
end

end