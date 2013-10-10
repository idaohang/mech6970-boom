function makemapped(h)
%MAKEMAPPED Convert ordinary graphics object to mapped object
%
% MAKEMAPPED(h) adds a Mapping Toolbox structure to the displayed objects 
% associated with h. h can be a handle, vector of handles, or any name string
% recognized by HANDLEM. The objects are then considered to be geographic 
% data. Objects extending outside the map frame should first be trimmed to the 
% map frame using TRIMCART.
%
% See also TRIMCART, HANDLEM, CART2GRN.

% needs something special for patches, so they can be reprojected

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.6.4.5 $  $Date: 2008/10/02 18:57:15 $

error(nargchk(1, 1, nargin, 'struct'))

if ischar(h) 
    h = handlem(h);
end

if ~ishghandle(h)
    error(['map:' mfilename ':notAGraphicsHandle'], ...
        'Inputs must be handles to graphics objects.')
end

% ensure vectors
h = h(:);

% Remove objects from the list that are already mapped
lengthin = length(h);
for i=length(h):-1:1
	if ismapped(h(i)); h(i) = []; end 
end

% Warn about them
if ~isequal(length(h), lengthin)
	warning('map:makemapped:objectAlreadyMapped', ...
        'Some objects already mapped.')
end

% Add a mapping toolbox object structure
set(h,'UserData',struct('trimmed',[],'clipped',[]),'ButtonDownFcn','uimaptbx');
