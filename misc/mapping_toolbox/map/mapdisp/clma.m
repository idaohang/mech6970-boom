function clma(action)
%CLMA Clear current map axes
%
%  CLMA clears the current map, deleting objects plotted on the map
%  but leaving the frame and grid lines displayed.
%
%  CLMA ALL clears the current map, frame and grid lines.  The map
%  definition is left in the axes definition.
%
%  CLMA PURGE removes the map definition from the current axes, but
%  leaves all objects projected on the axes intact.
%
%  See also HIDE, SHOWM, HANDLEM, NAMEM, TAGM, CLMO

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.9.4.3 $  $Date: 2007/11/09 20:26:55 $
% Written by:  E. Byrns, E. Brown

if nargin == 0
     action = 'map';
elseif ~strcmp(action,'map')   & ~strcmp(action,'all') & ~strcmp(action,'purge')
	validstr = strvcat('map','all','purge');
    indx     = strmatch(lower(action),validstr);

    if length(indx) == 1
       action = validstr(indx,:);
    else
       eid = sprintf('%s:%s:invalidCLMA', getcomp, mfilename);
       error(eid,'%s','Unrecognized CLMA action')
    end
end

%  Ensure that axes are for a map and get the current structure
mstruct = gcm;

%  Perform the appropriate deletion
switch action
      case 'map'
	          clmo('map')
			  mstruct.grid = 'off';
			  mstruct.parallellabel = 'off';
			  mstruct.meridianlabel = 'off';
              set(gca,'UserData',mstruct)
	  case 'all'
	          clmo('all')
			  mstruct.frame = 'off';
			  mstruct.grid  = 'off';
			  mstruct.parallellabel = 'off';
			  mstruct.meridianlabel = 'off';
              set(gca,'UserData',mstruct)
	  case 'purge'
	          showaxes on
			  set(gca,'UserData',[],...
			          'DataAspectRatioMode','auto',...
				 	  'ButtonDownFcn','');
end

