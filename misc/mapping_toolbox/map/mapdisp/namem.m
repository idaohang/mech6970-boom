function [obj,msg] = namem(hndl)
%NAMEM Determine names for valid map graphics objects
%
%  obj = NAMEM returns the object name for all objects on the current
%  axes.  The object name is defined as its tag, if the tag property
%  is supplied.  Otherwise, it is the object type.  Duplicate object
%  names are removed from the output string matrix.
%
%  obj = NAMEM(h) returns the object names for the objects specified
%  by the input handles h.
%
%  See also CLMO, HIDE, SHOWM, TAGM, HANDLEM.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.8.4.4 $  $Date: 2008/10/02 18:57:20 $

% Obsolete syntax
% ---------------
%  [obj,msg] = NAMEM(...)  returns a string msg indicating any error
%  encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Initialize inputs
if nargin == 0
    hndl = handlem('all');
end

%  Enforce a vector input
hndl = hndl(:);

%  Handle tests
if ischar(hndl) || any(~ishghandle(hndl))
    error(['map:' mfilename ':mapdispError'], ...
        'Vector of handles required');
end

%  Determine the types of each graphics object
%  Include in the name list only if it is unique
objtag  = get(hndl,'Tag');
objtype = get(hndl,'Type');

for i = 1:length(hndl)
    if length(hndl) == 1
        if isempty(objtag)
            objtag = objtype;
        end
    else
        if isempty(objtag{i})
            objtag{i} = objtype{i};
        end
    end
end

if isempty(hndl)
    obj = [];
elseif length(hndl) == 1
    obj = objtag;
else
    [obj,i] = unique(strvcat(objtag{:}),'rows');
    obj = strvcat(objtag{sort(i)});
end
