function schema
%SCHEMA Define the PolygonLayer class

%   Copyright 1996-2003 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2003/08/01 18:13:50 $

pkg = findpackage('MapModel');

% Extend the layer class
c = schema.class(pkg,'PolygonLayer',findclass(pkg,'layer'));

