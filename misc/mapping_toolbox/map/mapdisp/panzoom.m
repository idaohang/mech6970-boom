function panzoom(arg1, arg2, arg3)
%PANZOOM Pan and zoom on map axes
%
%  PANZOOM provides a graphical means for defining zoom limits
%  on a 2D plot.  A box which can be resized and moved is provided
%  on the current axes to allow the definition of a zoom area.  The
%  pan and zoom box can not be moved beyond the current axes limits
%
%  PANZOOM automatically toggles the current axes into a mode where
%  only actions recognized by PANZOOM are executed.  Upon exit of
%  this mode, all prior ButtonDown functions are restored to the
%  current axes and its children.
%
%  PANZOOM ON activates the pan and zoom tool.  PANZOOM OFF de-activates
%  the tool.  PANZOOM will toggle between these two states.
%
%  PANZOOM SETLIMITS will define the zoom out limits to the current
%  settings on the axes.  PANZOOM OUT will zoom out to the current axes
%  limit settings.  PANZOOM FULLVIEW will reset the axes to their
%  full view range and resets the pan and zoom tool with these settings.
%
%  PANZOOM recognizes the following keystrokes.  RETURN will set the axes
%  to the current zoom box and remain the in pan and zoom mode.  ENTER
%  sets the axes to the current zoom box and exits the pan and zoom
%  mode.  DELETE and ESCAPE will exit the pan and zoom mode (same as
%  PANZOOM OFF).
%
%  PANZOOM recognizes the following mouse actions when the cursor is inside
%  the pan and zoom box.  Single click and hold on center of box moves the
%  zoom box.  Single click and hold on corners allows the box to be resized.
%  Double click on the center of the box sets the axes to the current zoom
%  box and remains in pan and zoom mode (same as PANZOOM RETURN).
%
%  PANZOOM recognizes the following mouse actions when the cursor is
%  outside the pan and zoom box.  Single click relocates the pan and
%  zoom box centered on the cursor location.  Double click will
%  zoom out to the full view or last SETLIMITS setting.
%
%  Regardless of the cursor location, PANZOOM recognizes the following
%  mouse actions.  Extended click zooms out by a factor of 2.  Alternative
%  click exits the pan and zoom tool (same as PANZOOM OFF).
%
%  For Macintosh:   Extend click - Shift click mouse button
%                   Alternate click - Option click mouse button
%
%  For MS-Windows:  Extend click - Shift click left button or both buttons
%                   Alternate click - Control click left button or right button
%
%  For X-Windows:   Extend click - Shift click left button or middle button
%                   Alternate click - Control click left button or right button
%
%  See also ZOOM.

%  PANZOOM(AX, ...) is an undocumented syntax.  AX is the handle of the axes
%  object that PANZOOM should act on. AX may be followed by any of the
%  documented input syntaxes listed above.  This syntax was added to fix a bad
%  interaction between MAPTRIM and PANZOOM where excessive use of GCA in PANZOOM
%  caused a stray axes object to get created in a dialog window created by
%  MAPTRIM.  The syntax is undocumented because PANZOOM is in low-maintenance
%  mode and may be phased out eventually.  -SLE, 22-Oct-2009

% Copyright 1996-2009 The MathWorks, Inc.
% $Revision: 1.10.4.8 $  $Date: 2009/11/09 16:26:40 $
% Written by:  E. Byrns, E. Brown

%  Set default input to toggle mode

switch nargin
    case 0
        % panzoom()
        hndl = get(get(0, 'CurrentFigure'), 'CurrentAxes');
        action = 'toggle';
        
    case 1
        if ischar(arg1)
            % panzoom(action)
            action = arg1;
            hndl = get(get(0, 'CurrentFigure'), 'CurrentAxes');
        else
            % panzoom(ax)
            action = 'toggle';
            hndl = arg1;
        end
        
    case 2
        if ischar(arg1)
            % panzoom(action, buttonfcn)
            action = arg1;
            buttonfcn = arg2;
            hndl = get(get(0, 'CurrentFigure'), 'CurrentAxes');
        else
            % panzoom(ax, action)
            hndl = arg1;
            action = arg2;
        end
        
    case 3
        % panzoom(ax, action, buttonfcn)
        hndl = arg1;
        action = arg2;
        buttonfcn = arg3;
end    

% Ensure that the figure contains an axes
if isempty(hndl)
    error(['map:' mfilename ':noAxesInFigure'], 'No axes in current figure');
end

fig = ancestor(hndl, 'figure');

%  Ensure that the plot is in a 2D view

if any(get(hndl,'view') ~= [0 90]) && ~strcmp(action,'off')
    btn = questdlg(...
	       {'Must be in 2D view for operation.','Change to 2D view?'},...
		   'Incorrect View','Change','Cancel','Change');

    switch btn
	    case 'Change',      view(2);
		case 'Cancel',      return
    end
end


switch action
case 'on'        %  Display the pan&zoom box
       hbox = findobj(hndl,'Tag','panbox');   %  Delete old pan box if exists
       if ~isempty(hbox)                     %  But save any old function calls
	        fcns = get(hbox,'UserData');  delete(hbox)
	   else
	        fcns.windowdown = get(fig,'WindowButtonDownFcn');
	        fcns.windowup   = get(fig,'WindowButtonUpFcn');
	        fcns.windowmove = get(fig,'WindowButtonMotionFcn');
	        fcns.windowkey  = get(fig,'KeyPressFcn');
       end

       set(fig,'WindowButtonDownFcn','panzoom(''down'')',...  %  Set the
               'WindowButtonUpFcn','',...                     %  needed
               'WindowButtonMotionFcn','',...                 %  window
			   'KeyPressFcn','panzoom(''key'')')              %  functions


       xlim = get(hndl,'Xlim');            %  Compute the pan box coordinates
       xvec = min(xlim) + abs(diff(xlim)) * [1 2] / 3;
       xvec = xvec([1 2 2 1]);

       ylim = get(hndl,'Ylim');
	   yvec = min(ylim) + abs(diff(ylim)) * [1 2] / 3;
       yvec = yvec([2 2 1 1]);

	   zvec = max(get(hndl,'Zlim'));
	   zvec = zvec(ones(size(yvec)));

       boxcolor = 'red';       %  Pan box color

       ph = patch('Xdata',xvec,...    %  Display the pan box frame
           'Ydata',yvec, 'Zdata',zvec,...
           'Marker', 's',...
           'MarkerSize',3.5,...
           'MarkerFaceColor',boxcolor,...
           'MarkerEdgeColor',boxcolor,...
           'FaceColor','none',...
           'EdgeColor',boxcolor,...
           'Tag','panbox',...
           'UserData',fcns);
                   
       % Change state of zoom button in figure toolbar    
       hZmBtn = findall(fig,'TooltipString','Zoom In');
       if ~isempty(hZmBtn)
           set(hZmBtn,'State','on','ClickedCallback','panzoom(''zoomin'')');
       end
       % Change the font weight on the maptool zoom button if necessary
       hPshBtn = findobj(fig,'Type','uicontrol','Style','Push','String','Zoom');
       if ~isempty(hPshBtn)
           set(hPshBtn,'Fontweight','bold')
       end

       % re-set Window Button Motion Function
       axh = ancestor(ph, 'axes');
       set(fig,'WindowButtonMotionFcn',...
           @(varargin) panzoom(axh, 'zoomAffordances'));
       
case 'toggle'         %  Toggle pan&zoom box visible and invisible
       hbox = findobj(hndl,'Tag','panbox');
	   if isempty(hbox)           %  Initial call.  No box exists
	        panzoom('on')
	   elseif strcmp(get(hbox,'Visible'),'on')    %  Turn off
	        panzoom('off')
	   elseif strcmp(get(hbox,'Visible'),'off')    %  Turn on
	        panzoom('on')
	   end

case 'off'       %  Delete pan&zoom box and exit zoom mode
       hbox = findobj(hndl,'Tag','panbox');
       if ~isempty(hbox)
	       fcns = get(hbox,'UserData');  delete(hbox);
           set(fig,'WindowButtonDownFcn',fcns.windowdown,...   %  Reset all
	               'WindowButtonMotionFcn',fcns.windowmove,... %  window fcns
	               'WindowButtonUpFcn',fcns.windowup,...       %  which may have
			       'KeyPressFcn',fcns.windowkey)               %  been altered
       end
       % Change state of zoom button in figure toolbar    
       hZmBtn = findall(fig,'TooltipString','Zoom In');
       if ~isempty(hZmBtn)
           set(hZmBtn,'State','off','ClickedCallback','panzoom(''zoomin'')');
       end
       % Change the font weight on the maptool zoom button if necessary
       hPshBtn = findobj(fig,'Type','uicontrol','Style','Push','String','Zoom');
       if ~isempty(hPshBtn)
           set(hPshBtn,'Fontweight','normal')
       end

case 'zoomin'    %  Re-map ClickedCallback for Zoom Uitoogletool
    
       hZmBtn = findall(fig,'TooltipString','Zoom In');
       switch get(hZmBtn,'State')
       case 'on'
           zoom(fig,'on')
           % Change the font weight on the maptool zoom button if necessary
           hPshBtn = findobj(fig,'Type','uicontrol','Style','Push','String','Zoom');
           if ~isempty(hPshBtn)
               set(hPshBtn,'Fontweight','bold')
           end
       case 'off'
           zoom(fig,'off')
           % Change the font weight on the maptool zoom button if necessary
           hPshBtn = findobj(fig,'Type','uicontrol','Style','Push','String','Zoom');
           if ~isempty(hPshBtn)
               set(hPshBtn,'Fontweight','normal')
           end
           % Delete panbox if necessary
           hbox = findobj(hndl,'Tag','panbox');
           if ~isempty(hbox)
		       fcns = get(hbox,'UserData');  delete(hbox);
               set(fig,'WindowButtonDownFcn',fcns.windowdown,...   %  Reset all
		               'WindowButtonMotionFcn',fcns.windowmove,... %  window fcns
		               'WindowButtonUpFcn',fcns.windowup,...       %  which may have
				       'KeyPressFcn',fcns.windowkey)               %  been altered
           end
       end
       
case 'patch'     %  User has selected the pan&zoom box

        hbox   = findobj(hndl,'Tag','panbox');
		struct = get(hbox,'UserData');

        pt   = get(hndl,'CurrentPoint');              %  Determine the current
        xlim = get(hndl,'Xlim');                      %  point within the limits
		ylim = get(hndl,'Ylim');                      %  of the axes
        xpt  = max( min(xlim(2),pt(1,1)), xlim(1) );
        ypt  = max( min(ylim(2),pt(1,2)), ylim(1) );

        xvec = get(hbox,'Xdata');   % Remember that xvec and yvec are column
        yvec = get(hbox,'Ydata');   % vectors.  Necessary for patch movement below

		struct.xdel   = xvec - xpt;      %  Compute the delta of the current
		struct.ydel   = yvec - ypt;      %  point to the corners of the patch

		xeps = 0.15*abs(max(xvec) - min(xvec));   %  Corner hot range within
		yeps = 0.15*abs(max(yvec) - min(yvec));   %  15% of each corner


		if xeps <= 0.01*diff(xlim) || yeps <= 0.01*diff(ylim)
             xmiddle = min(xvec) + (max(xvec) - min(xvec))/2;    %  Box is too small
             ymiddle = min(yvec) + (max(yvec) - min(yvec))/2;    %  to grab a corner

             if        xpt <= xmiddle && ypt >  ymiddle;      struct.corner = 1;
                elseif xpt >  xmiddle && ypt >  ymiddle;      struct.corner = 2;
                elseif xpt >  xmiddle && ypt <= ymiddle;      struct.corner = 3;
                elseif xpt <= xmiddle && ypt <= ymiddle;      struct.corner = 4;
             end
		else
             cornertest = abs(struct.xdel) <= xeps & abs(struct.ydel) <= yeps;
             struct.corner = find(cornertest);     %  Identify which corner if any
		end

        if isempty(struct.corner)             %  Center of patch selected
		     struct.size     = get(hbox,'LineWidth');
             struct.sizeprop = 'LineWidth';
		     set(hbox,'LineWidth',2)
        else                                  %  Corner selected
		     struct.size     = get(hbox,'MarkerSize');
             struct.sizeprop = 'MarkerSize';
		     set(hbox,'MarkerSize',15)
        end

		set(fig,'WindowButtonMotionFcn','panzoom(''move'')',...  %  Movement
                'WindowButtonUpFcn','panzoom(''up'')')           %  callbacks

		set(hbox,'UserData',struct);   %  Save the needed data

case 'move'          %  Box is moving or resizing
        hbox   = findobj(hndl,'Tag','panbox');   %  Get box handle and
		struct = get(hbox,'UserData');          %  saved data

        pt   = get(hndl,'CurrentPoint');              %  Determine the current
        xlim = get(hndl,'Xlim');                      %  point within the limits
		ylim = get(hndl,'Ylim');                      %  of the axes
		xpt  = max( min(xlim(2),pt(1,1)), xlim(1) );
        ypt  = max( min(ylim(2),pt(1,2)), ylim(1) );

        xvec = get(hbox,'Xdata');   % Remember that xvec and yvec are column
        yvec = get(hbox,'Ydata');   % vectors.  Necessary for corner assignment below.

        if ~any(isnan(pt(:)))       %  Update only if cursor on axes.  PT is all NaNs if off axes
            if isempty(struct.corner)       %  Box movement
			     xvec = struct.xdel + xpt;     %  Compensate patch corners for
		         yvec = struct.ydel + ypt;     %  distance from cursor to edges

		    elseif struct.corner == 1           %  Corner resizing, upper left
                 xvec = [xpt; xvec([2;3]); xpt];    %  Make sure that xnew and ynew
			     yvec = [ypt([1;1]); yvec([3;4])];  %  are built as column vectors

		    elseif struct.corner == 2           %  Upper right
                 xvec = [xvec(1); xpt([1;1]); xvec(4)];
			     yvec = [ypt([1;1]); yvec([3;4])];

		    elseif struct.corner == 3           %  Lower right
                 xvec = [xvec(1); xpt([1;1]); xvec(4)];
			     yvec = [yvec([1;2]); ypt([1;1])];

		    elseif struct.corner == 4           %  Lower left
                 xvec = [xpt; xvec([2;3]); xpt];
			     yvec = [yvec([1;2]); ypt([1;1])];
            end
       end

	   set(hbox,'Xdata',xvec,'Ydata',yvec)   %  Update the box display

case 'up'        %  Mouse is released
       hbox   = findobj(hndl,'Tag','panbox');   %  Get box handle and
	   struct = get(hbox,'UserData');          %  saved data

       eval(['set(hbox,''',struct.sizeprop,''',struct.size)'])  %  Reset size

	   set(fig,'WindowButtonMotionFcn','panzoom(''zoomAffordances'')',... 
	           'WindowButtonUpFcn','')         %  Turn off movement callbacks

case 'apply'      %  Zoom to the current pan box

       hbox = findobj(hndl,'Tag','panbox');    %  Pan box and its data
       xdata = get(hbox,'Xdata');
	   ydata = get(hbox,'Ydata');

       set(hndl,'Xlim',[min(xdata) max(xdata)],...   %  Set axes limits
	             'Ylim',[min(ydata) max(ydata)])

	   panzoom('on');      %  Bring up a new pan&zoom box


case 'relocate'       %  Move the zoom box to the current mouse position

        hbox = findobj(hndl,'Type','patch','Tag','panbox');
        pt   = get(hndl,'CurrentPoint');              %  Determine the current
        xlim = get(hndl,'Xlim');                      %  point within the limits
		ylim = get(hndl,'Ylim');                      %  of the axes
        xpt  = max( min(xlim(2),pt(1,1)), xlim(1) );
        ypt  = max( min(ylim(2),pt(1,2)), ylim(1) );

        xvec = get(hbox,'Xdata');   % Remember that must xvec and yvec be column
        yvec = get(hbox,'Ydata');   % vectors.  Necessary for patch movement above

        delx = (max(xvec) - min(xvec))/2;      %  Compute 1/2 of each edge of
		dely = (max(yvec) - min(yvec))/2;      %  the pan box

		xlower = max(xlim(1),pt(1,1)-delx);    %  Compute the upper and lower box
		xupper = min(xlim(2),pt(1,1)+delx);    %  limits, centered about the
		ylower = max(ylim(1),pt(1,2)-dely);    %  current mouse point.  Ensure that
		yupper = min(ylim(2),pt(1,2)+dely);    %  the box lies within the axes

        xvec = [xupper; xlower; xlower; xupper];   %  Build the box coordinates
		yvec = [yupper; yupper; ylower; ylower];   %  as column vectors

        set(hbox,'Xdata',xvec,'Ydata',yvec)

case 'down'       %  Process the button down event properly

       hobject   = get(fig,'CurrentObject');        %  Current object
	   buttondwn = get(hobject,'ButtonDownFcn');    %  Potential callbacks
	   objtag    = get(hobject,'Tag');
	   selection = get(fig,'SelectionType');

       if strcmp(objtag,'panbox')        %  Click on pan box
          switch selection
			  case 'normal', panzoom('patch')
			  case 'open',   panzoom('apply')
			  case 'extend', zoom('down')        % Normal zoom operation
		      case 'alt',    panzoom('off')
           end

	   else
           switch selection
			  case 'normal',      panzoom('relocate')
			  case 'open',        zoom('out')
		      case 'extend',      zoom('down')
              case 'alt',         panzoom('off')
           end

%  In zoom mode, only the zoom button down actions (WindowButtonDown)
%  actions are to be carried out.  However, the button down event
%  that got you here will still be processed on the object's
%  ButtonDownFcn once the WindowButtonDown action is finished.
%  This is not the desired operation.  To get around this effect
%  the object ButtonDownFcn is redefined to a call back to
%  panzoom('buttondown').  This option of panzoom simply resets
%  the object ButtonDownFcn to its original status.  The original
%  status is identified by the second string argument to panzoom.
%  Thus, the original ButtonDown action of the individual object is
%  skipped, but not lost as the property of the object

              set(hobject,'ButtonDownFcn',...
	                      ['panzoom(''downreset'',''',buttondwn,''')'])

        end

case 'downreset'       %  Reset an object's button down function
      hobject = get(fig,'CurrentObject');
      set(hobject,'ButtonDownFcn',buttonfcn)


case 'key'             %  Key strokes recognized in zoom mode
	  keychar = abs(get(fig,'CurrentCharacter'));
	  if isempty(keychar); return; end
      switch keychar
	     case 3,      panzoom('apply');  panzoom('off')  %  Enter key
	     case 13,     panzoom('apply')                   %  Return key
	     case 27,     panzoom('off')                     %  Escape key
	     case 8,      panzoom('off')                     %  Delete key
	     case 127,    panzoom('off')                      %  Delete key (on PC)
     end

case 'out'       %  Zoom out to the current axes limits

	  zoom out            %  Clear the zoom limit settings
      panzoom('on')       %  Reinitialize the panzoom box

case 'fullview'       %  Reset an the display axes to the full view

      axis auto;          %  Reset the axes limits
	  zoom reset          %  Clear the zoom limit settings

      panzoom('on')       %  Reinitialize the panzoom box

case 'setlimits'          %  Set the zoom out limits to the current settings
	  zoom reset          %  Clear the zoom limit settings

case 'zoomAffordances'    % Pointer shape changes to indicate resizing or re-positioning
                          % of pan box.  A 'crosshair' pointer indicates that the pan
                          % box can be re-sized.  A 'fleur' pointer indicates that the
                          % pan box can be re-positioned.
   
	h = overobj2('axes');
    if ~isempty(h);
        hp = findobj(hndl,'Type','patch','Tag','panbox');
        fig = ancestor(hndl, 'figure');
        if ~isempty(hp)
            xd = get(hp,'XData'); yd = get(hp,'YData');
            mat = get(h,'CurrentPoint');
            cx = mat(1,1);
            cy = mat(1,2);
            xeps = 0.15*(max(xd)-min(xd));
            yeps = 0.15*(max(yd)-min(yd));
            xdel = abs(cx - xd);
            ydel = abs(cy - yd);
            if ~(min(xd) <= cx && cx <= max(xd) && min(yd) <= cy && cy <= max(yd))
                setptr(fig,'arrow');
            else
                if min(xdel) <= xeps && min(ydel) <= yeps
                    setptr(fig,'crosshair');
                else
                    setptr(fig,'fleur');
                end
            end
        end
    end
                          
end

% -------------------------------------------------------------------------------
% modified overobj to search for hidden objects also (line 478)

function h = overobj2(Type)

fig = get(0,'PointerWindow'); 
% Look for quick exit
if fig==0,
   h = [];
   return
end

% Assume root and figure units are pixels
p = get(0,'PointerLocation');
% Get figure position in pixels
figUnit = get(fig,'Units');
set(fig,'Units','pixels');
figPos = get(fig,'Position');
set(fig,'Units',figUnit)

x = (p(1)-figPos(1))/figPos(3);
y = (p(2)-figPos(2))/figPos(4);
% c = findobj(get(fig,'Children'),'flat','Type',Type,'Visible','on');
c = findobj(get(fig,'Children'),'flat','Type',Type);
for h = c',
   hUnit = get(h,'Units');
   set(h,'Units','norm')
   r = get(h,'Position');
   set(h,'Units',hUnit)
   if ( (x>r(1)) && (x<r(1)+r(3)) && (y>r(2)) && (y<r(2)+r(4)) )
      return
   end
end
h = [];
