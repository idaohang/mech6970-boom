function h = quiver3m(lat, lon, alt, u, v, w, varargin)
%QUIVER3M Project 3-D quiver plot on map axes
%
%  QUIVER3M(lat,lon,alt,u,v,w) projects three dimensional vector plot onto
%  the current map axes.  The vectors components (u,v,w) are specified at the
%  points (lat,lon,alt).  Note that both the (lat,lon) and (u,v) data must be
%  in the same angle units as the current map.  The units of alt and w must
%  be consistent, since these components will be used to define the altitude
%  of the vector above the map.  The vector is plotted from (lat,lon,alt) to
%  (lat+u,lon+v,alt+w).
%
%  QUIVER3M(lat,lon,alt,u,v,w,s) uses the input s to scale the velocity
%  vectors after they have been automatically scaled to fit within
%  the rectangular grid.  If omitted, s = 1 is assumed.  To suppress
%  the automatic scaling, specify s = 0.
%
%  QUIVER3M(lat,lon,alt,u,v,w,'LineSpec'),
%  QUIVER3M(lat,lon,alt,u,v,w,'LineSpec',s),
%  QUIVER3M(lat,lon,alt,u,v,w,'LineSpec','filled') and
%  QUIVER3M(lat,lon,alt,u,v,w,'LineSpec',s,'filled') use the LineSpec
%  string to define the line style of the velocity vectors.  If a symbol
%  is specified in 'LineSpec', then the symbol is plotted at the
%  base of the velocity vector.  Otherwise, an arrow is drawn at
%  the end of the velocity vector.  If a marker is specified at
%  the base, then this symbol can be filled in by providing
%  the input string 'filled'.
%
%  h = QUIVER3M(...) returns a vector of handles to the projected
%  velocity vectors.
%
%  See also  QUIVERM, QUIVER, QUIVER3.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.15.4.5 $  $Date: 2007/11/09 20:28:33 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
    errorOnMissingUI(mfilename)
end

%  Argument tests
if any([ndims(lat) ndims(lon) ndims(alt) ...
        ndims(u)   ndims(v)   ndims(w) ] > 2)
     error(['map:' mfilename ':inputContainsPages'], ...
         'Input data can not contain pages.')
     
elseif ~isequal(size(lat),size(lon),size(alt),size(u),size(v),size(w))
     error(['map:' mfilename ':inconsistentDims'], ...
         'Inconsistent dimensions for inputs.')
end

%  Parts of quiverm (quiver3m) closely parallel quiver (quiver3).
%  Unfortunately, you can not simply call quiver with projected
%  lat and lon data.  You do not get the appropriate clip and trim
%  data, which would preclude further re-projecting (using setm) of
%  the map.

%  Parse the optional input variables
switch length(varargin)
    case 1
        filled = [];
        if ~ischar(varargin{1});
            autoscale  = varargin{1};
            linespec = [];
        else
            linespec = varargin{1};
            autoscale = [];
        end

    case 2
        linespec = varargin{1};
        if ~ischar(varargin{2});
            autoscale  = varargin{2};
            filled = [];
        else
            filled = varargin{2};
            autoscale = [];
        end

    case 3
        linespec = varargin{1};
        autoscale = varargin{2};
        filled = varargin{3};

    otherwise
        linespec = [];
        autoscale = [];
        filled = [];
end

%  If unspecified, set autoscale to unity unless there
%  is only one unique point:
if isempty(autoscale)
    mstruct = gcm;
    if multipleDistinctLocations(lat,lon,mstruct.angleunits)
        autoscale = 1;
    else
        autoscale = 0;
    end
end

if ~isempty(linespec)
    [lstyle,lcolor,lmark,err] = colstyle(varargin{1});
    error(err) %#ok<ERTAG>
else
    lmark = [];
end

%  Autoscaling operation is taken directly from quiver3.
if autoscale
    % Base autoscale value on average spacing in the lat and lon
    % directions.  Estimate number of points in each direction as
    % either the size of the input arrays or the effective square
    % spacing if lat and lon are vectors.
    if min(size(lat))==1, n=sqrt(numel(lat));
        m=n;
    else
        [m,n]=size(lat);
    end
    delx = diff([min(lat(:)) max(lat(:))])/n;
    dely = diff([min(lon(:)) max(lon(:))])/m;
    delz = diff([min(alt(:)) max(alt(:))])/max(m,n);
    del  = sqrt(delx.^2 + dely.^2 + delz.^2);
    len  = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
    autoscale = autoscale*0.9 / max(len(:));
    u = u*autoscale;
    v = v*autoscale;
    w = w*autoscale;
end

%  Make inputs into row vectors.  Must be done after autoscaling
lat = lat(:)';
lon = lon(:)';
alt = alt(:)';
u   = u(:)';
v   = v(:)';
w   = w(:)';

%  Make the velocity vectors
vellat = [lat;  lat+u];
vellon = [lon;  lon+v];
velalt = [alt;  alt+w];
vellat(3,:) = NaN;
vellon(3,:) = NaN;
velalt(3,:) = NaN;

%  Set up for the next map
nextmap;

%  Project the velocity vectors as lines only
if ~isempty(linespec)
    h1 = linem(vellat(:),vellon(:),velalt(:),linespec,'Marker','none');
else
    h1 = linem(vellat(:),vellon(:),velalt(:),'Marker','none');
end

%  Make and plot the arrow heads if necessary
alpha = 0.33;   % Size of arrow head relative to the length of the vector
beta  = 0.33;   % Width of the base of the arrow head relative to the length
h2 = [];

if isempty(lmark)
    % Normalize beta
    beta = beta * sqrt(u.*u + v.*v + w.*w)/sqrt(u.*u + v.*v + eps*eps);

    % Make arrow heads and plot them
    hu = [lat+u-alpha*(u+beta*(v+eps));lat+u; ...
        lat+u-alpha*(u-beta*(v+eps))];       hu(4,:) = NaN;
    hv = [lon+v-alpha*(v-beta*(u+eps));lon+v; ...
        lon+v-alpha*(v+beta*(u+eps))];       hv(4,:) = NaN;
    hw = [alt+w-alpha*w;alt+w;alt+w-alpha*w];  hw(4,:) = NaN;

    if ~isempty(linespec)
        h2 = linem(hu(:),hv(:),hw(:),linespec,'Marker','none');
    else
        h2 = linem(hu(:),hv(:),hw(:),'Marker','none');
    end
end

%  Plot a marker on the base if necessary
h3 = [];
if ~isempty(lmark)
    h3 = linem(lat,lon,alt,linespec,'LineStyle','none');
    if strcmp(filled,'filled')
        set(h3,'MarkerFaceColor',get(h1,'color'));
    end
end

%  Set the tags
set([h1;h2;h3],'Tag','Quivers')

%  Set the output argument if necessary
if nargout == 1;
    h = [h1; h2; h3];
end
