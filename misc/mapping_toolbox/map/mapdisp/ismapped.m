function [tf,msg] = ismapped(h)
%ISMAPPED True if object is projected on map axes
%
%  TF = ISMAPPED returns 1 if the current object (gco) is projected on
%  a map axes.  Otherwise it returns 0.
%
%  TF = ISMAPPED(H) checks the object specified by the handle H.
%
%  [TF, MSG] = ISMAPPED(...) returns a string MSG explaining why the object
%  is not projected.
%
%  See also GCM, ISMAP.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.3 $  $Date: 2007/11/09 20:27:33 $

if nargin == 0
    h = get(get(0,'CurrentFigure'),'CurrentObject');
    if isempty(h)
        error(['map:' mfilename ':noObjectInFigure'], ...
            'No selected object in current figure.')
    end
end

%  Test the object parent for a map axes
parent = get(h,'Parent');
[isMapAxes,msg] = ismap(parent);
if isMapAxes
    % Validate object user data
    userdata = get(h,'UserData');
    if ~isstruct(userdata)
        tf = false;
        msg = 'Not a map object.';
    elseif ~all(isfield(userdata, {'trimmed','clipped'}))
        tf = false;
        msg = 'Not a map object.';
    else
        tf = true;
    end
else
    tf = false; 
end

% Preserve original behavior
tf = double(tf);
