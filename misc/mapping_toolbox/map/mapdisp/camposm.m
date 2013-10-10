function varargout = camposm(lat,long,alt)
%CAMPOSM Set camera position using geographic coordinates
%
%   CAMPOSM(lat,long,alt) sets the axes CameraPosition property of the
%   current map axes to the position specified in geographic coordinates. 
%   The inputs lat and long are assumed to be in the angle units of the 
%   current map axes. 
%
%   [x,y,z] = CAMPOSM(lat,long,alt) returns the camera position in the
%   projected cartesian coordinate system.
%
%   See also CAMTARGM, CAMUPM, CAMPOS, CAMVA

% Copyright 1996-2006 The MathWorks, Inc.
% $Revision: 1.7.4.2 $  $Date: 2006/06/15 20:11:57 $
% Written by: W. Stumpf, A. Kim, T. Debole

checknargin(3,3,nargin,mfilename);


[x,y,z] = mfwdtran(lat,long,alt);
set(gca,'CameraPosition',[x,y,z])

if nargout >= 1; varargout{1} = x; end
if nargout >= 2; varargout{2} = y; end
if nargout == 3; varargout{3} = z; end


