function [c, h] = contourfm(varargin)
%CONTOURFM  Project filled 2-D contour plot of map data
%
%   CONTOURFM(...) is the same as CONTOURM(...) except that the areas
%   between contours are filled with colors. For each contour interval,
%   CONTOURFM selects a distinct color from the figure's colormap.
%   You can obtain the same result by setting 'Fill','on' and
%   'LineColor','black' when calling CONTOURM.
%
%   Example 1
%   ---------
%   % Contour and fill the EGM96 geoid heights using 10 contour levels.
%   load geoid
%   worldmap world
%   contourfm(geoid, geoidrefvec, 10);
%
%   Example 2
%   ---------
%   % Contour and fill bathymetry and elevation data for the area around
%   % Korea, with contours at levels ranging from -5000 meters to 2500
%   % meters in increments of 500 meters.
%   load korea
%   figure
%   worldmap(map, refvec)
%   [c,h] = contourfm(map, refvec, -5000:500:2500);
%
%   See also CONTOURF, CONTOURM, CONTOUR3M.

% Copyright 1996-2010 The MathWorks, Inc.
% $Revision: 1.8.4.17 $  $Date: 2010/06/26 04:57:25 $

% Catch deprecated syntax.
if nargin == 0
    errorOnMissingUI('contourfm')
end

error(nargchk(2, inf, nargin, 'struct'))
switch(nargout)
    case 0
        contourm(varargin{:},'Fill','on','DefaultLineColor','black');
    case 1
        c = contourm(varargin{:},'Fill','on','DefaultLineColor','black');
    case 2
        [c, h] = contourm(varargin{:},'Fill','on','DefaultLineColor','black');
end

% Note: The 'DefaultLineColor' parameter that appears above is for
% internal use only and may be subject to change. Outside the toolbox,
% the 'LineColor' parameter should be used with contourm instead.
