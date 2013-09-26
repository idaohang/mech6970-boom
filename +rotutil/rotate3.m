function R = rotate3(varargin)
%% ROTATE3 creates rotation matrices given the 3-axis rotation angles (eulerian)
% Rotate first about the x-axis, then the y-axis, and finally about the
% z-axis.
% 
% USAGE:
%   R = rotate3(x,y,z, units)
%   R = rotate3([x,y,z], units)
% 
% INPUT:
%   units: either 'deg' or 'rad'
%   x,y,z: double
%
% TODO:
%   - length checking for inputs - do 3-D arrays of output

if nargin == 2
    xyz = varargin{1}; x = xyz(1); y = xyz(2); z = xyz(3);
    units = varargin{2};
elseif nargin == 4
    x = varargin{1}; y = varargin{2}; z = varargin{3};
    units = varargin{4};
end

if strcmp(units,'deg')
    x = rad2deg(x); y = rad2deg(y); z = rad2deg(z); 
elseif strcmp(units,'rad')
    %
else
    error('need "rad" or "deg" for units');
end
    
Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];

R = Rz*Ry*Rx;
    
% end

end