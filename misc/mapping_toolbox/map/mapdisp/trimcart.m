function trimcart(h)
%TRIMCART Trim graphic objects to map frame
%
%  TRIMCART(h) trims the cartesian graphic objects to the map frame.  h can 
%  be a handle or a vector of handles to graphics objects.  h can also be any 
%  object name recognized by HANDLEM. TRIMCART trims lines, surfaces, and 
%  text objects.  TRIMCART does not trim patches.
%  
%  See also HANDLEM, MAKEMAPPED.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.5.4.5 $  $Date: 2008/10/02 18:57:36 $
% Written by:  W. Stumpf

error(nargchk(1, 1, nargin, 'struct'))
if ischar(h)
    h = handlem(h);
end

if ~ishghandle(h)
    error(['map:' mfilename ':notAHandle'], ...
        'Inputs must be handles to graphics objects');
end

% inpolygon seems optimized for speed over memory efficiency. Reduce the 
% amount of data it requires by thinning the frame.

hframe = handlem('Frame');
newframe = 0;
if isempty(hframe)
	hframe = framem;
	newframe = 1;
end

xframe = get(hframe,'Xdata');
yframe = get(hframe,'Ydata');   

[xframe,yframe] = reducem(xframe(1:end-1),yframe(1:end-1)); % 0,0 at end of frame

% Loop over the objects

for i=1:length(h)
		switch get(h(i),'type')
	
		case 'line'
			clipcartline(h(i),xframe,yframe)
		case 'surface'
			clipcartsurf(h(i),xframe,yframe)
		case 'text'
			clipcarttext(h(i),xframe,yframe)
		otherwise
			warning(['map:' mfilename ':cannotClipToFrame'], ...
                'Can''t clip %s to the map frame', get(h(i),'type'))
		end
end

% to toggle frame state back to what it was

if newframe; delete(hframe); end  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clipcartline(h,xframe,yframe)

% data for the line object

xline = get(h,'xdata');
yline = get(h,'ydata');
zline = get(h,'zdata');
if isempty(zline); zline = zeros(size(xline)); end

% clip
in = +inpolygon(xline,yline,xframe,yframe);
in(~in) = NaN;
in(~isnan(in)) = 1;

xline = xline.*in; 
yline = yline.*in;
zline = zline.*in;

% restore

set(h,'xdata',xline,'ydata',yline,'zdata',zline);

			
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clipcartsurf(h,xframe,yframe)


% data for the surface object

xgrid = get(h,'xdata');
ygrid = get(h,'ydata');

% ensure grid is a matrix

if min(size(xgrid)) && min(size(ygrid))
	[xgrid,ygrid] = meshgrid(xgrid,ygrid);
end

% inpolygon is optimized for speed over memory usage. 
% Take grid in chunks to avoid running out of memory.

for i=1:size(ygrid,1)
	if any(~isnan(ygrid(i,:)))
		in = +inpolygon(xgrid(i,:),ygrid(i,:),xframe,yframe);
		in(~in) = NaN;
		in(~isnan(in)) = 1;
		xgrid(i,:) = xgrid(i,:).*in; ygrid(i,:) = ygrid(i,:).*in;
	end
end

set(h,'xdata',xgrid,'ydata',ygrid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clipcarttext(h,xframe,yframe)


% data for the line object

pos = get(h,'Position');

x = pos(1);
y = pos(2);

% clip
in = inpolygon(x,y,xframe,yframe);

if ~in
	set(h,'Visible','off');
end

