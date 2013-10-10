function [tf,msg] = ismap(varargin)
%ISMAP  True for axes with map projection
%
%  TF = ISMAP returns 1 if the current axes (gca) is a map axes.
%  Otherwise, it returns 0.
%
%  TF = ISMAP(H) checks the axes specified by the handle H.
%
%  [TF, MSG] = ISMAP(...) returns a string MSG explaining why a map axes
%  was not found. 
%
%  See also GCM, ISMAPPED.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.10.4.4 $  $Date: 2007/11/09 20:27:31 $

try
    gcm(varargin{:});
    tf = true;
    msg = '';
catch e
    if strcmp(e.identifier, 'map:gcm:noAxesInFigure')
        throw(e)
    end
    tf = false;
    msg = e.message;
end
    
% Preserve original behavior
tf = double(tf);
