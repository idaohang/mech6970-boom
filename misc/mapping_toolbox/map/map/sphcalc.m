function val = sphcalc(rad,calculation)

%SPHCALC  Computes volume and surface area for a sphere
%
%  SPHCALC(r,'volume') computes the volume of a sphere
%  defined by the input radius.  The units are defined
%  by the input radius.
%
%  SPHCALC(r,'surfarea') computes the surface area of a sphere
%  defined by the input radius.
%
%  See also ELPCALC.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.2 $  $Date: 2007/11/09 20:25:34 $
% Written by:  E. Byrns, E. Brown

if max(size(rad)) ~= 1
   error(['map:' mfilename ':mapError'], 'Radius input must be a scalar')
end

if strcmp(calculation,'volume')
    val = (4*pi/3) * rad^3;

elseif strcmp(calculation,'surfarea')
    val = 4*pi*rad^2;

else
    error(['map:' mfilename ':mapError'], ...
        'Unrecognized calculation string')
end
