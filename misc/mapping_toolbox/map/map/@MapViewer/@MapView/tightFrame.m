function tightFrame(this)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.2 $  $Date: 2008/10/26 14:26:41 $

% The toolbar and display panel need about 70 pixels, so the height needs
% 60 added to it.

[width, height] = this.Axis.getVisibleSizeInPixels();

figurePosition = getpixelposition(this.Figure);
figurePosition(3) = max(width, this.MinWidth);
figurePosition(4) = max(height,this.MinHeight) + 70;
setpixelposition(this.Figure,figurePosition)
