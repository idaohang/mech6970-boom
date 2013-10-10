function h = RasterLayer(name)
%RASTERLAYER constructor for raster image layer object.
%
%   RASTERLAYER(NAME) constructs a raster layer with the name NAME.  There
%   are no components initially in the raster layer.  Components can be
%   added with the addcomponent method.  RGBComponents or
%   IntensityComponents may be added to a RasterLayer object. However, the
%   layer can only contain one type of component. 

%   Copyright 1996-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.4 $  $Date: 2006/12/10 20:04:29 $

h = MapModel.RasterLayer;

h.ComponentType = {'RGBComponent','IntensityComponent',...
                   'IndexedComponent','GriddedComponent'};
h.LayerName = name;
h.Visible = 'on';