function previewmap
%PREVIEWMAP Preview map at printed size
% 
%   PREVIEWMAP changes the size of the current figure to match the printed 
%   output.  This provides an accurate display of the relative placement and 
%   size of objects as they will be printed.  If the resulting figure size 
%   exceeds the screen size, the figure will be enlarged as much as possible.
%   
%   See also PAPERSCALE, AXESSCALE.

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.4.4.4 $  $Date: 2008/05/14 22:01:27 $

% Written by: W. Stumpf

h = gcf;

figureunits = get(h,'units');
paperunits = get(h,'paperunits');paperpos = get(h,'paperpos');
set(h,'units',paperunits,'pos',paperpos)

pos = get(h,'Position');

if max(abs(pos(3:4) - paperpos(3:4))) > 0.05
	set(h,'units',paperunits,'pos',[0 0 paperpos(3:4)])
	if max(abs(pos(3:4) - paperpos(3:4))) > 0.05
		warning('map:previewmap:figureExceedsScreenSize',...
            'Figure window larger than screen size. Enlarging as much as possible')
	end
end

set(h,'units',figureunits)

shiftwin(gcf)

figure(gcf)
