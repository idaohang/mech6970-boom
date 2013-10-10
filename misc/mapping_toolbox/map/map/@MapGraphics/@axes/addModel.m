function addModel(this,model)

% Copyright 1996-2008 The MathWorks, Inc.
% $Revision: 1.1.6.8 $  $Date: 2009/11/09 16:25:48 $

this.Model = model;

% Install Listeners
this.Listeners = [this.Listeners(:);...
                  handle.listener(model,'LayerAdded',{@newlayer model this});...
                  handle.listener(model,'LayerRemoved',{@removelayer this});...
                  handle.listener(model,'LayerOrderChanged',{@layerorderchanged this});...
                  handle.listener(model,'ShowBoundingBox',{@showBoundingBox this});...
                  handle.listener(model,'Visible',{@setVisible this})
                 ];

% Render the graphics in the model
model.render(this.getAxes());

%-------------------------------------------------------------------------
function h = createComponentListener(layer,model,hMapAxes)
componentProp = findprop(layer,'Components');
h = handle.listener(layer,componentProp,...
                    'PropertyPostSet',{@newcomponent model hMapAxes});

%-------------------------------------------------------------------------%
function reorderChildren(hMapAxes,layerorder)
newChildren = [];
for i=1:length(layerorder)
  newChildren = [newChildren; hMapAxes.getLayerHandles([layerorder{i} '_BoundingBox']);...
                 hMapAxes.getLayerHandles(layerorder{i})];
end
assert(numel(newChildren) == numel(get(hMapAxes.getAxes(),'Children')), ...
    'MapGraphics:axes:addModel:layerHandleMismatch',...
    'Number of re-ordered children fails to match original count.')
set(hMapAxes.getAxes(),'Children',newChildren);
refresh(ancestor(hMapAxes.getAxes(),'Figure'))

%%%%%%%%%%%%%%%%%%%LISTENERS%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
function setVisible(src,eventData,hMapAxes) %#ok<INUSL>
set(hMapAxes.getLayerHandles(eventData.Name),...
    'Visible',eventData.Value);

%-------------------------------------------------------------------------%
function showBoundingBox(src,eventData,hMapAxes) %#ok<INUSL>
delete(hMapAxes.getLayerHandles([eventData.Name '_BoundingBox']));
hMapAxes.Model.getLayer(eventData.Name).renderBoundingBox(hMapAxes.getAxes());


%-------------------------------------------------------------------------%
function layerorderchanged(src,eventData,hMapAxes) %#ok<INUSL>
reorderChildren(hMapAxes,eventData.layerorder);

%-------------------------------------------------------------------------%
function removelayer(src,eventData,hMapAxes) %#ok<INUSL>
% Remove Layer
delete(hMapAxes.getLayerHandles(eventData.LayerName));
% Remove Bounding Box
delete(hMapAxes.getLayerHandles([eventData.LayerName '_BoundingBox']))

%-------------------------------------------------------------------------%
function newlayer(src,eventData,model,hMapAxes) %#ok<INUSL>
% Update graphics when a layer is added or removed
layername = eventData.LayerName;
model.renderLayer(hMapAxes.getAxes(),layername);
% Add a listener to the new layer's component property
hMapAxes.ComponentListener= createComponentListener(model.getLayer(layername),model,hMapAxes);

%-------------------------------------------------------------------------%
function newcomponent(src,eventData,model,hMapAxes) %#ok<INUSL>
% Update graphics when a new component is added to a layer
layer = eventData.AffectedObject;
component = eventData.NewValue(end,:);
layer.renderComponent(hMapAxes.getAxes(),component);
reorderChildren(hMapAxes,model.getLayerOrder);

