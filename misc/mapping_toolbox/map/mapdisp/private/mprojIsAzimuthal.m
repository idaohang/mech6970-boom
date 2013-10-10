function tf = mprojIsAzimuthal(mapprojection)
% True if the projection with the ID string corresponding to
% MAPPROJECTION is azimuthal or pseudo-azimuthal.

% Copyright 2009 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2009/05/14 17:06:37 $

tf = any(strcmp(mapprojection, {...
       'ups','stereo','ortho','breusing',...
       'eqaazim','eqdazim','gnomonic','vperspec','wiechel'}));
