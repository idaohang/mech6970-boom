function mstruct = initmstruct
% Create an empty mstruct with minimal initialization.

% Copyright 2008 The MathWorks, Inc.
% $Revision: 1.1.6.1 $  $Date: 2008/05/12 21:33:14 $

defaultColor = get(0,'DefaultAxesXcolor');
mstruct = struct( ...
    'mapprojection', [], ...
    'zone',          [], ...
    'angleunits',    'degrees', ...
    'aspect',        'normal', ...
    'falsenorthing', [], ...
    'falseeasting',  [], ...
    'fixedorient',   [], ...
    'geoid',         [1 0], ...
    'maplatlimit',   [], ...
    'maplonlimit',   [], ...
    'mapparallels',  [], ...
    'nparallels',    0, ...
    ...
    'origin',        [], ...
    'scalefactor',   [], ...
    'trimlat',       [], ...
    'trimlon',       [], ...
    ...
    'frame',         [], ...
    'ffill',         100, ...
    'fedgecolor',    defaultColor, ...
    'ffacecolor',    'none', ...
    'flatlimit',     [], ...
    'flinewidth',    2, ...
    'flonlimit',     [], ...
    ...
    'grid',          [], ...
    'galtitude',     inf, ... % Signal to put graticule at upper end of zlim
    'gcolor',        defaultColor, ...
    'glinestyle',    ':', ...
    'glinewidth',    0.5, ...
    ...
    'mlineexception', [], ...
    'mlinefill',      100, ...
    'mlinelimit',     [], ...
    'mlinelocation',  [], ...
    'mlinevisible',   'on', ...
    ...
    'plineexception', [], ...
    'plinefill',      100, ...
    'plinelimit',     [], ...
    'plinelocation',  [], ...
    'plinevisible',   'on', ...
    ...
    'fontangle',      'normal', ...
    'fontcolor',      defaultColor, ...
    'fontname',       get(0,'factoryAxesFontName'), ...
    'fontsize',       get(0,'factoryAxesFontSize'), ...
    'fontunits',      get(0,'factoryAxesFontUnits'), ...
    'fontweight',     get(0,'factoryAxesFontWeight'), ...
    'labelformat',    'compass', ...
    'labelrotation',  'off', ...
    'labelunits',     [], ...
    ...
    'meridianlabel',  [], ...
    'mlabellocation', [], ...
    'mlabelparallel', [], ...
    'mlabelround',    0, ...
    ...
    'parallellabel',  [], ...
    'plabellocation', [], ...
    'plabelmeridian', [], ...
    'plabelround',    0);
