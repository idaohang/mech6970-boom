function setVisible(this,b)
%SETVISIBLE
%
%   SETVISIBLE(B) Sets the graphics components to be visible (B = 'On' or true)
%   or invisible (B = 'Off' or false).

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.3 $  $Date: 2008/11/24 14:59:05 $

if islogical(b)
  if b
    b = 'on';
  else
    b = 'off';
  end
end

set(this.LineHandle,'Visible',b)
set(this.TextHandle,'Visible',b)
