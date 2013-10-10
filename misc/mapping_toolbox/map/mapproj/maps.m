function [idstr,msg] = maps(projin)
%MAPS  List available map projections and verify names
%
%   MAPS lists the available map projections.
%
%   s = MAPS IDLIST returns the list of all available map
%   projection id strings.
%
%   s = MAPS NAMELIST returns the list of all available map
%   projection names.
%
%   s = MAPS CLASSCODES returns the list of codes for the projection
%   classes of all available maps projections.
%
%   str = MAPS('str') verifies and standardizes the map projection
%   id string.
%
%   See also AXESM.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.1.6.4 $  $Date: 2007/11/09 20:30:03 $
% Written by:  E. Byrns, E. Brown

% Obsolete Syntax
% ---------------
%   [str,msg] = MAPS(...) returns the string indicating any error
%   condition encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Get the map list structure
%  The map list structure is of the form:
%        list.Name
%        list.IdString
%        list.Classification
%        list.ClassCode

list = maplist;

%  Initialize outputs if necessary

if nargout ~= 0
    idstr = [];
end

if nargin == 0
    %  Display available projections
    formatstr = '%-20s  %-32s    %-15s  \n';
    fprintf('\n%-20s \n\n','MapTools Projections')
    fprintf(formatstr,'CLASS','NAME','ID STRING');
    for i = 1:length(list)
        fprintf(formatstr,list(i).Classification,...
            list(i).Name,...
            list(i).IdString);
    end
    fprintf('\n%s\n','* Denotes availability for sphere only')
    fprintf('\n\n');
elseif strcmp(projin,'namelist')
    %  Return the list of available map names
    idstr = strvcat(list(:).Name);
    indx = find(idstr == '*');    % Eliminate the asterisks in the names
    if ~isempty(indx)
        idstr(indx) = ' ';
    end
elseif strcmp(projin,'idlist')
    %  Return the list of available id strings
    idstr = strvcat(list(:).IdString);
elseif strcmp(projin,'classcodes')
    %  Return the list of available id strings
    idstr = strvcat(list(:).ClassCode);
elseif ischar(projin)
    %  Convert from string name to projection number
    projin = projin(:)';     %  Enforce row string vector
    idstring = {list(:).IdString};
    idx = strmatch(projin,idstring);  %  Test for a match
    if numel(idx) == 1
        idstr = idstring{idx};   %  Set the name string
    elseif numel(idx) > 1
        % Remove the non-matching elements from idstring
        idstring = idstring(idx);
        k = find(strcmpi(projin,idstring));
        if numel(k) == 1
            idstr = idstring{k};
        else
            error(['map:' mfilename ':mapprojError'], ...
                'Nonunique projection name:  %s',projin)
        end
    else
        error(['map:' mfilename ':mapprojError'], ...
            'Projection name not found:  %s',projin)
    end
else
    %  Not a valid option
    error(['map:' mfilename ':mapprojError'], ...
        'Unrecognized projection string')
end
