function mat = getm(hndl,propname)
%GETM Map object properties
%
%  mat = GETM(h,'MapPropertyName') returns the value of the specified
%  map property for the map graphics object with handle h.  The
%  graphics object h must be a map axis or one of its children.
%
%  mat = GETM(h) returns all map property values for the map object with
%  handle h.  
%
%  GETM MAPPROJECTION lists the available map projections.
%  GETM AXES lists the map axes properties.
%  GETM UNITS lists the recognized unit strings.
%
%  See also AXESM, SETM, GET, SET.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.11.4.4 $  $Date: 2008/10/02 18:57:11 $
% Written by:  E. Byrns, E. Brown

%  Programmers Note:  GETM is a time consuming call because of
%  the strmatch and getfield functions.  It is advised that when
%  programming, direct access is made to the map structure and
%  do not use calls to GETM.  For an example of this, see GRIDM
%  and FRAMEM.

error(nargchk(1, 2, nargin, 'struct'))

if nargin == 1
    %  Handle provided.  Get all properties
    if ~ischar(hndl)
	     propname = [];
    else
        %  No handle.  Test for special string
        validparm = {'mapprojection','units','axes'};
		indx = strmatch(lower(hndl),validparm);
        if length(indx) == 1
		     hndl = validparm{indx};
		else
		    error(['map:' mfilename ':invalidString'], ...
                'Unrecognized GETM string.')
        end

         switch  hndl
		    case 'mapprojection'
                maps;
                return
		    case 'units'
                unitstr;
                return
		    case 'axes'
                setm(gca); 
                return
         end
    end
end

%  Valid handle test
if isempty(hndl)
    error(['map:' mfilename ':emptyHandle'], 'Handle is not valid.')
end
if ~isscalar(hndl)
    error(['map:' mfilename ':multipleHandles'], ...
        'Multiple handles not allowed.')
end
if ~ishghandle(hndl)
    error(['map:' mfilename ':invalidHandle'], 'Handle is not valid.')
end

%  Get the corresponding axis handle
if ishghandle(hndl,'axes')
     maphndl = hndl;
elseif ishghandle(get(hndl,'Parent'),'axes')
     maphndl = get(hndl,'Parent');
else
     error(['map:' mfilename ':handleNotAxis'], ...
         'GETM only works with a Map Axis and its children.')
end

%  Test for a valid map axes and get the corresponding map structure
gcm(maphndl);

%  Get the user data structure from this object
userstruct = get(hndl,'UserData');
if ~isstruct(userstruct)
    error(['map:' mfilename ':expectedStruct'], ...
        'Map structure not found in object')
end

%  Return the entire structure if propname is empty
if isempty(propname)
    mat = userstruct;
    return
end

%  Otherwise, get the fields of the structure and test for a match
structfields = fieldnames(userstruct);
indx = strmatch(lower(propname),lower(structfields));
if isempty(indx)
    error(['map:' mfilename ':invalidProperty'], ...
        'Incorrect property for object.')
elseif length(indx) == 1
    propname = structfields{indx};
else
	indx = strmatch(lower(propname),lower(structfields),'exact');	  
	if length(indx) == 1
    	propname = structfields{indx};
	else
	    error(['map:' mfilename ':nonUniqueProperty'], ...
       'Property %s name not unique - supply more characters.', propname)
	end
end

%  If match is found, then return the corresponding property
mat = userstruct.(propname);
