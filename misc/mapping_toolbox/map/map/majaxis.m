function semimajor = majaxis(in1,in2)
%MAJAXIS  Semimajor axis of ellipse with given semiminor axis and
%         eccentricity
%
%  a = MAJAXIS(semiminor,e) computes the semimajor axis of an ellipse
%  (or ellipsoid of revolution) given the semiminor axis and eccentricity.
%  The input data can be scalar or matrices of equal dimensions.
%
%  a = MAJAXIS(vec) assumes a 2 element vector (vec) is supplied,
%  where vec = [semiminor, e].
%
%  See also AXES2ECC, FLAT2ECC, MINAXIS, N2ECC.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.4 $  $Date: 2007/11/09 20:24:28 $
% Written by:  E. Byrns, E. Brown

error(nargchk(1, 2, nargin, 'struct'))

if nargin == 1
    if ~isequal(sort(size(in1)),[1 2])
        error(['map:' mfilename ':mapError'], ...
            'Input must be a 2 element vector')
    else
        in1 = ignoreComplex(in1, mfilename, 'vec');
        semiminor = in1(1);
        eccent    = in1(2);
    end
elseif nargin == 2
    if ~isequal(size(in1),size(in2))
        error(['map:' mfilename ':mapError'], ...
            'Inconsistent input dimensions')
    else
        semiminor = ignoreComplex(in1, mfilename, 'semiminor');
        eccent    = ignoreComplex(in2, mfilename, 'eccentricity');
    end
end

%  Compute the semimajor axis
semimajor = semiminor ./ sqrt(1 - eccent.^2);
