function [lat,lon,alt] = cart2grn(hndl,mstruct)
%CART2GRN Transform projected coordinates to Greenwich system
%
%  [lat,lon,alt] = CART2GRN will transform the current displayed
%  map object from cartesian coordinates to the Greenwich frame.  In
%  doing so, any clips or trims introduced during the display process
%  are removed from the output data.
%
%  [lat,lon,alt] = CART2GRN(hndl) will transform the object specified
%  by hndl, which must currently be displayed on a map axes.
%
%  [lat,lon,alt] = CART2GRN(hndl,mstruct) uses the input map
%  structure to determine the projection parameters for the
%  displayed object.  If this is omitted, then the structure from
%  the current axes will be obtained.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.8.4.3 $  $Date: 2008/10/02 18:57:04 $
% Written by:  E. Byrns, E. Brown


%  Argument tests

if nargin == 0
    hndl = gco;
    mstruct = [];
elseif nargin == 1
    mstruct = [];
end

%  Ensure that a single object of a map axes if provided

if ~isscalar(hndl) || ~ishghandle(hndl)
     eid = sprintf('%s:%s:invalidObjectHandle', getcomp, mfilename);
     error(eid,'Scalar object handle required.');
elseif ~ismapped(hndl)
     eid = sprintf('%s:%s:invalidObject', getcomp, mfilename);
     error(eid,'Object must be mapped.')
elseif ~ismap(get(hndl,'Parent'))
    eid = sprintf('%s:%s:noObjectParent', getcomp, mfilename);
    error(eid,'Object must be on a map axes.')
end

%  Get the map structure if necessary

if isempty(mstruct)
    mstruct = gcm(get(hndl,'Parent'));
end

%  Get the object's map coordinates and height

if ishghandle(hndl,'text')
    
    pos = get(hndl,'Position');
    x = pos(1);
    y = pos(2);
    if length(pos) == 2
        z = 0;
    else
        z = pos(3);
    end
    
elseif ishghandle(hndl,'light')
    
    pos = get(hndl,'Position');
    x = pos(1);
    y = pos(2);
    if length(pos) == 2
        z = 0;
    else
        z = pos(3);
    end
    
elseif ishghandle(hndl,'patch')
    
    vertices = get(hndl,'Vertices');
    x = vertices(:,1);
    y = vertices(:,2);
    if size(vertices,2) == 3;
        z = vertices(:,3);
    else
        z = zeros(size(x));
    end
    
else   %  Just get the x,y and z data for all other objects
    
    x = get(hndl,'Xdata');
    y = get(hndl,'Ydata');
    z = get(hndl,'Zdata');
    
    if size(x,1) == 1            %  Ensure column vectors
        x = x';
        y = y';
        z = z';
    end
end

%  Get the clip and trim data

savepts = get(hndl,'UserData');

%  Transform from the projected space to Greenwich space

[lat,lon,alt] = minvtran(mstruct,x,y,z,get(hndl,'Type'),savepts);
