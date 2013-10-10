function ecc = axes2ecc(in1,in2)
%AXES2ECC  Eccentricity of ellipse with given axes lengths
%
%   ECC = AXES2ECC(SEMIMAJOR,SEMIMINOR) computes the eccentricity
%   of an ellipse (or ellipsoid of revolution) given the semimajor
%   and semiminor axes.  The input data can be scalar or matrices
%   of equal dimensions.
%
%   ECC = AXES2ECC(VEC) assumes a 2 element vector (VEC) is supplied,
%   where VEC = [SEMIMAJOR SEMIMINOR].
%
%   See also  MAJAXIS, MINAXIS, ECC2FLAT, ECC2N.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.9.4.4 $    $Date: 2010/09/24 14:33:08 $

error(nargchk(1,2,nargin,'struct'))

if nargin == 1
    if numel(in1) ~= 2
        error('map:axes2ecc:numelVecNot2','VEC must be a 2 element vector.')
    else
        validateattributes(in1,{'double'},{'real','positive'},mfilename,'VEC',1);
        semimajor = in1(1);
        semiminor = in1(2);
    end
elseif nargin == 2
    semimajor = in1;
    semiminor = in2;
    validateattributes(semimajor,{'double'},{'real','positive'},mfilename,'SEMIMAJOR',1);
    validateattributes(semiminor,{'double'},{'real','positive'},mfilename,'SEMIMINOR',2);
    if ~isequal(size(semimajor),size(semiminor))
        error('map:axes2ecc:inconsistentSizes', 'Inconsistent input sizes.')
    end
end

if any(semiminor > semimajor)
    error('map:axes2ecc:invalidAxes', ...
        'SEMIMINOR axis may not be larger than SEMIMAJOR axis.');
end

%  Compute the eccentricity
ecc = sqrt(semimajor.^2 - semiminor.^2) ./ semimajor;
