function [hout,htout] = mdistort(action,levels,gsize)
% MDISTORT Display contours of constant map distortion
% 
%   MDISTORT, with no input arguments, toggles the display of contours
%   of projection-induced distortion on the current map axes.  The
%   magnitude of the distortion is reported in percent.
%   
%   MDISTORT OFF removes the contours.
%   
%   MDISTORT('PARAMETER') or MDISTORT PARAMETER displays contours of 
%   distortion for the specified parameter.  Recognized parameter
%   strings are 'area', 'angles' for the maximum angular distortion of
%   right angles, 'scale' or 'maxscale' for the maximum scale,
%   'minscale' for the minimum scale, 'parscale' for scale along the
%   parallels, 'merscale' for scale along the meridians, and
%   'scaleratio' for the ratio of maximum and minimum scale.  If
%   omitted, the 'maxscale' parameter is displayed.  All parameters are
%   displayed as percent distortion, except angles, which are displayed
%   in degrees.
% 
%   MDISTORT('PARAMETER',LEVELS) specifies the levels for which the
%   contours are drawn.  LEVELS is a vector of values as used by
%   CONTOUR. If omitted, the default levels are used.
%   
%   MDISTORT('PARAMETER',LEVELS,GSIZE) controls the size of the
%   underlying graticule used to compute the contours.  GSIZE is a
%   two-element vector containing the number of rows and columns.  If
%   omitted, the default Mapping Toolbox graticule size of [50 100] is
%   assumed.
% 
%   h = MDISTORT(...) returns a handle to the contourgroup object
%   containing the contours and text.
%
%   Note:  MDISTORT may not be valid when used with UTM, and distortion
%          is minimal within a given UTM zone. MDISTORT issues a warning
%          if a UTM projection is encountered.
%
%   See also TISSOT, DISTORTCALC, VFWDTRAN.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.5.4.12 $  $Date: 2010/03/22 03:52:12 $
% Written by: W. Stumpf

% Reference
% ---------
% Maling, Coordinate Systems and Map Projections, 2nd Edition, 

error(nargchk(0, 3, nargin, 'struct'))
ax = gca;

if nargin < 1
	h = findobj(ax,'Tag','DISTORTMlines');
    if isempty(h)
		action = 'maxscale';
	else
		action = 'off';
    end
end

if nargin < 2
    levels = [];
end

if nargin < 3
	gsize = [];
end

if strcmpi(action,'off')
    % Return the name of the currently plotted parameter for redisplay
    % after a projection change
    h = mdistortOff(ax);
    ht = [];
else
    [h, ht] = mdistortConstruct(ax, action,levels,gsize);
end

% Set output arguments
if nargout >= 1;
    hout = h;
end

if nargout == 2
    htout = ht;
    warning('map:mdistort:textHandlesRequested', ...
        ['Text object handles are now available as children of', ...
        ' the contourgroup returned by %s, and use of the', ...
        ' optional second output %s is no longer recommended.'], ...
        'MDISTORT', 'ht')
end

%-----------------------------------------------------------------------

function [h, ht] = mdistortConstruct(ax, action, levels, gsize)

% Check that action is a string
if ~ischar(action)
    error(['map:' mfilename ':argNotString'], ...
        'Argument must be a string.')
end
action = lower(action);

% Check for recognized actions
actions = {'on','off','area','angles','angle','scale','maxscale',...
           'minscale','parscale','merscale','scaleratio'};

if isempty(strmatch(action,actions,'exact'))
	error(['map:' mfilename ':unknownOption'], 'Unknown MDISTORT option.')
end

switch action
case 'on'
	action = 'maxscale';
end

% Get the map structure from the axes
mstruct = getm(ax);

% Issue a warning if projection is UTM
if strmatch(mstruct.mapprojection,'utm')
    warning('map:mdistort:projectionIsUTM',...
        'MDISTORT is not intended for use with UTM');
end

% Construct a graticule within the frame limits and compute 
% distortion parameters at the graticule intersections.

% Grid up the area within the frame and project to map coordinates
[xgrat, ygrat] = framegrat(mstruct, gsize);

% Convert from planar (x-y) system to geographic latitude-longitude
[latgrat, longrat] = minvtran(mstruct, xgrat, ygrat);

% Assign scaling parameter
param = setParameter(action, mstruct, latgrat, longrat);

% Set up contour levels
if isempty(levels)
    if any(strmatch(action, {'angles','angle'}, 'exact'))
        levels = [ 0 0.5 1 2 5 10 20 30 45 90 180];
    else
        levels = [ 0 0.5 1 2 5 10 15 20 50 100  200 400 800];
        levels = [-fliplr(levels) levels(2:end)]; % used for scale calculations
    end
end

% Remove any previously plotted results
mdistortOff(ax);

% Add contour lines and labels to the plot
ht = [];
[c,h] = contour(ax, xgrat, ygrat, param, levels);
if ~isempty(h); 
	tagm(h,'DISTORTMlines')
	ht = clabel(c,h);

	if ~isempty(ht)
		set(ht,'color','r')
		tagm(ht,'DISTORTMtext')
	end
	
    s = struct('clipped', [], 'trimmed', [], 'parameter', action);
	set([h;ht],'UserData',s)
end

%-----------------------------------------------------------------------

function param = mdistortOff(ax)

% Note: The text objects (with handles in array ht above) are children
% of h, so there's no need to delete them explicitly.
h = findobj(ax,'Tag','DISTORTMlines');
param = [];
if ~isempty(h)
    struct = getm(h(1));
    param = struct.parameter;
end
delete(h)

%-----------------------------------------------------------------------

function param = setParameter(action, mstruct, latgrat, longrat)

% Compute the projection distortion parameters by finite differences
[areascale, angdef, maxscale, minscale, merscale, parscale] = ...
	distortcalc(mstruct, latgrat, longrat);

switch action
    case 'area'
        param = (areascale-1)*100; % in percent
    case {'angles','angle'}
        % Convert angular deformation to degrees
        param = toDegrees(mstruct.angleunits, angdef);
    case {'maxscale','scale'}
        param = (maxscale-1)*100; % in percent
    case 'minscale'
        param = (minscale-1)*100; % in percent
    case 'parscale'
        param = (parscale-1)*100; % in percent
    case 'merscale'
        param = (merscale-1)*100; % in percent
    case 'scaleratio'
        param = (abs(minscale./maxscale)-1)*100; % in percent
end

%-----------------------------------------------------------------------

function [xgrat, ygrat] = framegrat(mstruct, gsize)
% Construct a graticule that covers the map frame in the frame's own
% latitude-longitude system (which may be shifted and rotated with
% respect to geographic coordinates), and project it to map X-Y.

%  Reset the projection origin, moving to the frame's system.
projImplementedViaRotation = ...
    ~any(strcmp(mstruct.mapprojection, {'tranmerc', 'cassinistd', ...
            'eqaconicstd', 'eqdconicstd', 'lambertstd', 'polyconstd'}));
        
if projImplementedViaRotation
    mstruct.origin = [0 0 0];
else
    % Don't modify the origin latitude -- this projection does not
    % simply rotate an auxiliary sphere.
    mstruct.origin = [mstruct.origin(1) 0 0];
end

% Construct the graticule.
epsilon = 100000*epsm('degrees');
if ~mprojIsAzimuthal(mstruct.mapprojection)
    % non-azimuthal frame
    flatlim = mstruct.flatlimit + epsilon * [1 -1];
    flonlim = mstruct.flonlimit + epsilon * [1 -1];
	[latgrat, longrat] = meshgrat(flatlim, flonlim, gsize);
else
    % azimuthal frame
	rnglim = mstruct.flatlimit;
	rnglim(1) = 0;
	azlim = fromDegrees(mstruct.angleunits, [0+epsilon 360-epsilon]);	
	[rnggrat,azgrat] = meshgrat(rnglim,azlim, gsize);
	[latgrat, longrat] = reckon( ...
        'gc', 0, 0, rnggrat, azgrat, mstruct.angleunits);
end

% Project the graticule.
[xgrat, ygrat] = mfwdtran(mstruct, latgrat, longrat);
