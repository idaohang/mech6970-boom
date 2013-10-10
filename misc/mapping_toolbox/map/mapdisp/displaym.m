function [h,msg] = displaym(varargin)
%DISPLAYM Display geographic data from display structure
%
%  DISPLAYM(DISPLAY_STRUCT) projects and displays the data contained in
%  a Mapping Toolbox display structure, DISPLAY_STRUCT, in the current
%  axes.  The current axes must be a map axes.
%
%  DISPLAYM(DISPLAY_STRUCT, STR) displays the vector data elements of
%  DISPLAY_STRUCT whose 'tag' fields contains strings beginning with the
%  string STR.  Vector data elements are those whose 'type' field is
%  either 'line' or 'patch'.  The string match is case-insensitive.
%
%  DISPLAYM(DISPLAY_STRUCT, STRINGS) displays the vector data elements
%  of DISPLAY_STRUCT whose 'tag' field matches begins with one of the
%  elements (or rows) of STRINGS.  STRINGS is a cell array of strings
%  (or a 2-D character array).  In the case of character array, trailing
%  blanks are stripped from each row before matching.
%
%  DISPLAYM(DISPLAY_STRUCT, STRINGS, SEARCHMETHOD) controls the method
%  used to match the values of the 'tag' field in DISPLAY_STRUCT, as
%  follows:
%
%     'strmatch'   Search for matches at the beginning of the tag
%                  (similar to the STRMATCH function)
%
%     'findstr'    Search within the tag (similar to the FINDSTR function)
%
%     'exact'      Search for exact matches
%
%  Note that when SEARCHMETHOD is specified the search is case-sensitive.
%
%  h = DISPLAYM(DISPLAY_STRUCT) returns handles to the graphic objects
%  created by DISPLAYM.
%
%  Note
%  ----
%  Use GEOSHOW instead of DISPLAYM if you need to display a geostruct
%  (created using SHAPEREAD or GSHHS, for example).
%
% See also EXTRACTM, GEOSHOW, MLAYERS, UPDATEGEOSTRUCT.

% Copyright 1996-2007 The MathWorks, Inc.
% $Revision: 1.15.4.5 $  $Date: 2007/12/10 21:39:13 $
% Written by:  E. Byrns, E. Brown, W. Stumpf

% Obsolete syntax
% ---------------
% [h,msg] = DISPLAYM(...) returns a string indicating any error encountered.
if nargout > 1
    warnObsoleteMSGSyntax(mfilename)
    msg = '';
end

%  Initialize output if necessary
if nargout ~= 0
    h = [];
end

% handle cell arrays from rootlayr or vmap0ui
if nargin >= 1 && iscell(varargin{1}) && size(varargin{1},2) == 2
    h0 = [];
    c=varargin{1};
    for i=1:size(c,1)
        if nargin == 1
            hi = displaym(c{i,1});
        else
            hi = displaym(c{i,1},varargin{2:end});
        end

        if ~isempty(hi)
            h0 = [h0 hi];
        end
    end
    return
end

error(nargchk(1, 3, nargin, 'struct'))

if nargin == 1;
    mstruct = varargin{1};
    if isempty(mstruct)
        return
    end
elseif nargin == 2 || nargin == 3;
    mstruct = varargin{1};varargin(1) = [];
    if isempty(mstruct)
        return
    end
    warnstate = warning('off','map:extractm:ignoringNonvectorData');
    [lat,lon,indx] = extractm(mstruct,varargin{:});
    warning(warnstate)
    h0 = displaym(mstruct(indx));
    if nargout > 0
        h = h0;
    end
    return
end

colororder = get(gca,'ColorOrder');          %  lines if no properties supplied
cdataoffset = length(get(gca,'Children'));   %  Indexing of patch face color

utags = unique(strvcat(mstruct.tag),'rows'); %  unique tags, for coloring

%  Loop through the structure, plotting each member
h0 = [];

for i = 1:length(mstruct)
    if isempty(mstruct(i).altitude);
        switch mstruct(i).type
            case {'patch','regular'}
                alt = 0;
            otherwise
                alt = 0*ones(size(mstruct(i).lat));
        end

    else
        alt = mstruct(i).altitude;
    end

    %  Plot using other properties if they're supplied.
    %  If applying properties causes an error, try again without. This may
    %  result in some residue on the plot.

    if ~isempty(mstruct(i).otherproperty) && iscell(mstruct(i).otherproperty)
        switch mstruct(i).type
            case 'light'
                try
                    htemp = lightm(mstruct(i).lat,mstruct(i).long,alt,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = lightm(mstruct(i).lat,mstruct(i).long,alt);
                end

            case 'line'

                try
                    htemp = linem(mstruct(i).lat,mstruct(i).long,alt,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = linem(mstruct(i).lat,mstruct(i).long,alt);
                end
            case 'patch'
                try
                    htemp = patchesm(mstruct(i).lat,mstruct(i).long,alt,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = patchesm(mstruct(i).lat,mstruct(i).long,alt);
                end
            case 'regular'
                try
                    htemp = meshm(mstruct(i).map,mstruct(i).maplegend,...
                        mstruct(i).meshgrat,alt,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = meshm(mstruct(i).map,mstruct(i).maplegend,...
                        mstruct(i).meshgrat,alt);
                end
            case 'surface'
                try
                    htemp = surfacem(mstruct(i).lat,mstruct(i).long,...
                        mstruct(i).map,alt,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = surfacem(mstruct(i).lat,mstruct(i).long,...
                        mstruct(i).map,alt);
                end
            case 'text'
                try
                    htemp = textm(mstruct(i).lat,mstruct(i).long,alt,...
                        mstruct(i).string,...
                        mstruct(i).otherproperty{:});
                catch
                    htemp = textm(mstruct(i).lat,mstruct(i).long,alt,...
                        mstruct(i).string);
                end
        end
    else     %  No other properties
        switch mstruct(i).type
            case 'light'
                htemp = lightm(mstruct(i).lat,mstruct(i).long,alt);
            case 'line'
                % to color all similarly tagged objects the same
                indx = strmatch(mstruct(i).tag,utags,'exact');
                % no match, tag was empty
                if isempty(indx); indx = i;  end
                %  Adjust the default color spec
                clrcount = 1+mod(indx,size(colororder,1)) ;
                clrstring = colororder(clrcount,:);

                htemp = linem(mstruct(i).lat,mstruct(i).long,alt,...
                    'Color',clrstring);
            case 'patch'     %  Use face color indexing as a default

                % to color all similarly tagged objects the same
                indx = strmatch(mstruct(i).tag,utags,'exact');
                if isempty(indx); indx = i;  end            % no match, tag was empty

                htemp = patchesm(mstruct(i).lat,mstruct(i).long,alt,...
                    'Cdata',cdataoffset+indx,'FaceColor','flat');
            case 'regular'
                htemp = meshm(mstruct(i).map,mstruct(i).maplegend,...
                    mstruct(i).meshgrat,alt);
            case 'surface'
                htemp = surfacem(mstruct(i).lat,mstruct(i).long,...
                    mstruct(i).map,alt);
            case 'text'
                htemp = textm(mstruct(i).lat,mstruct(i).long, ...
                    alt,...
                    mstruct(i).string);
        end
    end

    h0 = [h0 htemp(:)'];
    set(htemp,'Tag',mstruct(i).tag)
end

%  Set output arguments if necessary
if nargout > 0
    h = h0;
end
