function varargs = importFromFileAndSetDataArgs(fcnName, varargin)
%importFromFileAndSetDataArgs Read data from file if filename provided
%
%   If the first element of VARARGIN is a string, assume that it's a data
%   file name, import the data, insert the data at the beginning of
%   VARARGIN and add a DisplayType parameter/value pair to correspond to
%   the type of data file; otherwise return the argument list unaltered.
%
%   fcnName is just a string containing the function name to be used only
%   in creating error messages.
%
%   See also GEOSHOW, MAPSHOW.

% Copyright 2006 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2006/06/15 20:11:04 $

firstArgIsFilename = ~isempty(varargin) && ischar(varargin{1});
if firstArgIsFilename 
   filename = varargin{1};
   varargin(1) = [];
   
   % Import the data from the file into the cell array dataArgs.
   % For example, when reading from a SDTS raster DEM file, dataArgs will
   % contain two elements.  The first element is a 2-D array containing an
   % elevation grid and the second element is the referencing matrix.
   % defaultDisplayType is a string containing the default display type for
   % this type of data file.
   [dataArgs, defaultDisplayType] = readMapData(fcnName, filename);
   
   % Place the default display type near the beginning of the argument list
   % so that it may be overridden from the command line (expressed as a
   % subsequent parameter/value pair in varargin).
   varargs = {dataArgs{:}, 'DisplayType', defaultDisplayType, varargin{:}};
else
   varargs = varargin;
end
