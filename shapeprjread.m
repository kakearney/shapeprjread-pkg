function [Shp, m, punitfac] = shapeprjread(file, varargin)
%SHAPEPRJREAD Read a projected shapefile
%
% Shp = shapeprjread(file)
% Shp = shapeprjread(file, p1, v1, ...)
% [Shp, m] = shapeprjread(file, p1, v1, ...)
%
% This function reads data from a shapefile, incorporating information from
% a projection file.  If the projection file indicates a geographic
% coordinate system, data will be read into a geographic data structure.
% If the projection file indicates a projected coordinate system, the
% original coordinates (X,Y) will be read and then reverse projected to
% Lat/Lon coordinates.
%
% Input variables:
%
%   file:   name of file, either the base name or any of the component
%           files accepted by shaperead (.shp, .shx, or .dbf files).  A
%           corresponding .prj file should be located in the same location.
%
% Optional input variables (passed as parameter/value pairs)
%
%   see shaperead.m.  All parameters except 'UseGeoCoords' will be
%   accepted.  Note that all selectors, including BoundingBox, will be
%   applied to the data as it is stored in the file, not to the
%   reverse-projected data.
%
% Output variables:
%
%   Shp:    structure holding combination mapstruct/geostruct (includes
%           both X/Y and Lat/Lon fields, plus any other attribute fields
%           found in the file)
%
%   m:      map projection structure used for reverse projection.
%
%   fac:    projection unit conversion factor.  When call minvtran with the
%           above map projection structure, multiple x/y values by this
%           first.  When going the opposite way, divide after calling
%           mfwdtran.

%-------------------------
% Setup
%-------------------------

% Parse filename

[pth, base, ext] = fileparts(file);
prjfile = fullfile(pth, [base '.prj']);

%-------------------------
% Read data from files
%-------------------------

if ~exist(prjfile)
    warning('No projection file found; reading as mapstruct');
    Shp = shaperead(file, varargin{:});
else
    Prj = prjread(prjfile);
    
    % Determine which coordinate system is used
    
    if isfield(Prj, 'COMPD_CS') % Compound coordinate system
        
        if isfield(Prj.COMPD_CS, 'PROJCS')
           isproj = true;
            P = Prj.COMPD_CS.PROJCS;
        else
            isproj = false;
        end
        
    elseif isfield(Prj, 'PROJCS') % Projected coordinate system
        
        P = Prj.PROJCS;
        isproj = true;
        
    elseif isfield(Prj, 'GEOCCS') % Geocentric coordinate system
        
        error('Geocentric coordinates not yet supported');
        
    elseif isfield(Prj, 'GEOGCS') % Geographic coordinate system
        
        isproj = false;
        
    end
    
    % Read shapefile data, and if necessary, reverse project
    
    if isproj
        [m, punitfac] = buildmstruct(P);
        
        Shp = shaperead(file, varargin{:});
        for ii = 1:length(Shp)
            [Shp(ii).Lat, Shp(ii).Lon] = minvtran(m, Shp(ii).X*punitfac, Shp(ii).Y*punitfac);
        end
    else
        Shp = shaperead(file, 'UseGeoCoords', true, varargin{:});
    end
    
end    

%-------------------------
% Create map projection
% structure
%-------------------------

function [m, punitfac] = buildmstruct(P)

% Projections: ESRI name and Matlab name (I can't find this documented, so
% I'll just have to build it up manually)

projections = {...
    'Transverse_Mercator'           'tranmerc'
    'Albers'                        'eqaconicstd'
    'Albers Conical Equal Area'     'eqaconicstd'
    'Lambert_Azimuthal_Equal_Area'  'eqaazim'
};

[tf,loc] = ismember(P.PROJECTION.name, projections(:,1));
if ~tf
    error('Need to add this projection to the list: %s', P.PROJECTION.name);
end 

% Convert into mstruct

if loc == 1
    if regexpfound(P.name, 'UTM') % No idea if this is robust or not
        m = defaultm('utm');
        zone = regexp(P.name, 'UTM_Zone_([A-Z|0-9]*)', 'tokens', 'once');
        zone = zone{1};
        m.zone = zone;
    else
        m = defaultm(projections{loc,2});
    end
else
    m = defaultm(projections{loc,2});
end

keywords = {...
'AUTHORITY'    
'AXIS'         
'COMPD_CS'     
'CONCAT_MT'    
'DATUM'        
'FITTED_CS'    
'GEOCCS'       
'GEOGCS'       
'INVERSE_MT'   
'LOCAL_DATUM'  
'LOCAL_CS'     
'PARAMETER'    
'PARAM_MT'
'PASSTHROUGH_MT'
'PRIMEM'      
'PROJCS'      
'PROJECTION'  
'SPHEROID'    
'TOWGS84'     
'UNIT'        
'VERT_DATUM'  
'VERT_CS'};

flds = fieldnames(P);
isparam = ~ismember(flds, [keywords; 'name']);
flds = flds(isparam);

% Get units first, since it will scale other properties

if isfield(P, 'UNIT') % This applies only to linear distances in parameters
    punitfac = P.UNIT.conversionFactor;
else
    punitfac = 1;
end

if isfield(P.GEOGCS, 'UNIT') % This applies to angular distances in geoid
    gunitfac = P.GEOGCS.UNIT.conversionFactor;
else
    gunitfac = 1;
end

% Set the geoid

m.geoid = oblateSpheroid;
m.geoid.SemimajorAxis = P.GEOGCS.DATUM.SPHEROID.semiMajorAxis; % Assume in m?
m.geoid.InverseFlattening = P.GEOGCS.DATUM.SPHEROID.inverseFlattening;

% Set remaining parameters if they are found

for ii = 1:length(flds)
    switch flds{ii}
        case {'False_Easting', 'false_easting'}
            m.falseeasting = str2num(P.(flds{ii})) * punitfac;
        case {'False_Northing', 'false_northing'}
            m.falsenorthing = str2num(P.(flds{ii})) * punitfac;
        case {'Central_Meridian', 'central_meridian'}
            m.origin(2) = str2num(P.(flds{ii}));
        case {'Latitude_Of_Origin', 'Central_Parallel', 'latitude_of_origin'}
            m.origin(1) = str2num(P.(flds{ii}));
        case 'Scale_Factor'
            m.scale = str2num(P.(flds{ii}));
        case 'Standard_Parallel'
            m.mapparallel = str2num(P.(flds{ii}));
        case {'Standard_Parallel_1', 'standard_parallel_1'}
            m.mapparallels(1) = str2num(P.(flds{ii}));
        case {'Standard_Parallel_2', 'standard_parallel_2'}
            m.mapparallels(2) = str2num(P.(flds{ii}));
        otherwise
            warning('Found extra parameter: %s', flds{ii}); % May need to update code for these
    end
end

m = defaultm(m);




