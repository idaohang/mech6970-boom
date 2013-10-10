function restack(hobj,action)
%RESTACK Restack objects within map axes
%
%   RESTACK(h,position) changes the stacking order of the object h
%   within the axes.  h can be a handle or vector of handles to graphics
%   objects, or h can be a name string recognized by HANDLEM. Recognized
%   position strings are 'top','bottom','bot','up' or 'down'. RESTACK
%   permutes the order of the children of the axes.
%
%   See also MOBJECTS, ZDATAM.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.6.4.4 $  $Date: 2008/10/02 18:57:26 $

% Check for recognized actions
actions = {'top','bottom','bot','up','down'};
if isempty(strmatch(action,actions,'exact'))
	error(['map:' mfilename ':unknownPosition'], ...
        'Unknown %s option.', 'RESTACK')
end

% Get the handle of the target object
if ~ishghandle(hobj)
	hobj = handlem(hobj);
end

if isempty(hobj)
    return
end

% Get the handle of the map axes
ax = ancestor(hobj(1),'axes');

% Get the handles of the children of the axes
children = get(ax,'Children');

% identify the location of the target object within the children
[ignored,objpos] = intersect(children,hobj);
objpos = sort(objpos);

% identify the locations everything besides the target object
[ignored,notobjpos] = setxor(children,hobj);
notobjpos = sort(notobjpos);

% reorder the children 
switch action
case 'top'	
	newchildren = [children(objpos(:)); children(notobjpos(:))];
    
case {'bottom','bot'}
	newchildren = [children(notobjpos(:)); children(objpos(:))];
    
case 'down'	
	topobj = max(objpos);
	newchildren = [	children( notobjpos( notobjpos(:)<= (topobj+1) ) ); ...
					children(objpos(:)); ...
					children( notobjpos( notobjpos(:)>topobj+1 ) ) ];

case 'up'	
	botobj = min(objpos);
	newchildren = [	children( notobjpos( notobjpos(:)< (botobj-1) ) ); ...
					children(objpos(:)); ...
					children( notobjpos( notobjpos(:)>= (botobj-1) ) ) ];

end

set(ax,'Children',newchildren)
