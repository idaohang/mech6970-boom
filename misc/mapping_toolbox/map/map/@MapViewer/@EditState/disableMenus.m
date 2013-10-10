function disableMenus(this)

%   Copyright 1996-2005 The MathWorks, Inc.
%   $Revision: 1.1.6.3 $  $Date: 2005/03/31 16:33:42 $

% Turn off all edit menu items
set(get(this.EditMenu,'Children'),'Enable','off');

