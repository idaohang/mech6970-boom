function axesscale(baseaxes,ax)
%AXESSCALE Resize axes for equivalent scale
% 
%   AXESSCALE resizes all axes in the current figure to have the same scale 
%   as the current axes (gca). In this context, scale means the
%   relationship between axes X and Y coordinates to figure and paper
%   coordinates. The XLimMode and YLimMode of the axes are set to Manual to
%   prevent autoscaling from changing the scale
% 
%   AXESSCALE(hbase) uses the axes hbase as the reference axes, and
%   rescales the other axes in the current figure.
% 
%   AXESSCALE(hbase,hother) uses the axes hbase as the base axes, and
%   rescales only the axes hother.
%
%   Example
%   -------
%   Display the conterminous United States, Alaska, and Hawaii in separate
%   axes in the same figure, with a common scale.
%
%   % Read state names and coordinates, extract Alaska and Hawaii
%   states = shaperead('usastatehi', 'UseGeoCoords', true);
%   statenames = {states.Name};
%   alaska = states(strmatch('Alaska', statenames));
%   hawaii = states(strmatch('Hawaii', statenames));
%
%   % Create a figure for the conterminous states
%   f1 = figure; hconus = usamap('conus'); tightmap
%   geoshow(states, 'FaceColor', [0.5 1 0.5]);
%   load conus gtlakelat gtlakelon
%   geoshow(gtlakelat, gtlakelon,...
%           'DisplayType', 'polygon', 'FaceColor', 'cyan')
%   framem off; gridm off; mlabel off; plabel off
%
%   % Working figure for additional calls to usamap
%   f2 = figure('Visible','off');
%
%   halaska = axes; usamap('alaska'); tightmap;
%   geoshow(alaska, 'FaceColor', [0.5 1 0.5]);
%   framem off; gridm off; mlabel off; plabel off
%   set(halaska,'Parent',f1)
%
%   hhawaii = axes; usamap('hawaii'); tightmap; 
%   geoshow(hawaii, 'FaceColor', [0.5 1 0.5]);
%   framem off; gridm off; mlabel off; plabel off
%   set(hhawaii,'Parent',f1)
%
%   close(f2)
% 
%   % Arrange the axes as desired
%   set(hconus,'Position',[0.1 0.25 0.85 0.6])
%   set(halaska,'Position',[0.019531 -0.020833 0.2 0.2])
%   set(hhawaii,'Position',[0.5 0 .2 .2])
%
%   % Resize alaska and hawaii axes
%   axesscale(hconus)
%   hidem([halaska hhawaii])
% 
%   See also PAPERSCALE.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.6.4.6 $  $Date: 2008/10/02 18:57:03 $
% Written by: W. Stumpf, L. Job

checknargin(0,2,nargin,mfilename);
if nargin == 0
	baseaxes = gca;
	ax = findobj(gcf,'type','axes');
elseif nargin == 1
	ax = findobj(gcf,'type','axes');
end

% check that the base handle is to an axes
if length(baseaxes(:)) > 1
   eid = sprintf('%s:%s:invalidBaseAxesLength', getcomp, mfilename);
   error(eid,'Only one base axes at a time.'); 
end

if any(~ishghandle([baseaxes;ax(:)])) 
   eid = sprintf('%s:%s:invalidHandle', getcomp, mfilename);
   error(eid,'Inputs must be axis handles.');
end

if ~ishghandle(baseaxes,'axes')
   eid = sprintf('%s:%s:invalidBaseAxesType', getcomp, mfilename);
   error(eid,'Inputs must be axis handles.');
end

ellipsoid = []; 
if ismap(baseaxes)
	ellipsoid = getm(baseaxes,'geoid');
end

% check that the other handles are to axes

warned = 0;
for i=1:length(ax)
    if ~ishghandle(ax(i),'axes')
        eid = sprintf('%s:%s:invalidHandle', getcomp, mfilename);
        error(eid,'Inputs must be axis handles.');
    end

% Check that ellipsoids are approximately the same to ensure that scaling 
% between geographic and x-y data is the same

	if ismap(ax(i))
		thisellipsoid = getm(baseaxes,'geoid');
	else
		thisellipsoid = [];
	end
	
	if ~isempty(ellipsoid) && ~isempty(thisellipsoid) && abs((ellipsoid(1)-thisellipsoid(1))/ellipsoid(1))  > 0.01 ; 
        if ~warned
            wid = sprintf('%s:%s:ellipsoidUnitsDiffer', getcomp, mfilename);
            warning(wid, ...
                ['Ellipsoids differ between axes.' ...
                ' Check that units of ellipsoids are the same.'])
            warned = 1;
        end
	end
end

% get the properties of BASEAXES

xlim = get(baseaxes,'xlim');
ylim = get(baseaxes,'ylim');
pos = get(baseaxes,'pos');
deltax = pos(3);
deltay = pos(4);

% set the xscale and yscale

xscale = deltax/abs(diff(xlim));
yscale = deltay/abs(diff(ylim));

% loop on remaining axes

for i = 1:size(ax,1)
    if ishghandle(ax(i),'axes')
        xlim = abs(diff(get(ax(i),'xlim')));
        ylim = abs(diff(get(ax(i),'ylim')));
        pos = get(ax(i),'pos');
        pos(3) = xscale*xlim;
        pos(4) = yscale*ylim;
        set(ax(i),'pos',pos);
    end
end

% Lock down xlim and ylim to ensure that scale remains constant if
% additional data is plotted

set([baseaxes;ax(:)],'xLimMode','manual','ylimmode','manual');
