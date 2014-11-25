function S = prjread(file)
%PRJREAD Read ESRI shapefile projection file
%
% S = prjread(file)
%
% Read in the projection component of a shapefile (which Matlab's shaperead
% ignores).
%
% Input variables:
%
%   file:   .prj file
%
% Output files:
%
%   S:      structure, often nested, holding all data from the file.  Field
%           names correspond to either keywords (hold structures) or
%           properties (hold values)



% Copright 2012 Kelly Kearney


fid = fopen(file, 'rt');
txt = textscan(fid, '%s', 'delimiter', '\n');
txt = cat(2, txt{1}{:});
fclose(fid);


S = parsewkt(txt);

% Subfunction: Parse the well-known text format into a structure

function S = parsewkt(str)

A = parsegroup(str);

switch A.key
    
    case 'AUTHORITY'    % ["<name>", "<code>"]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).code = A.elements{2};
 
    case 'AXIS'         % ["<name>", NORTH | SOUTH | EAST | WEST | UP | DOWN | OTHER]
        
        name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        name = ['AXIS_' name];
        S.(name) = A.elements{2};
        
%         S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
%         S.(A.key).orient = A.elements{2};
        
    case 'COMPD_CS'     % ["<name>", <head cs>, <tail cs> {,<authority>}]
    
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
    
    case 'CONCAT_MT'
    case 'DATUM'        % ["<name>", <spheroid> {,<to wgs84>} {,<authority>}]
    
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
    
    case 'FITTED_CS'    % ["<name>", <to base>, <base cs>]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:3
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'GEOCCS'       % ["<name>", <datum>, <prime meridian>, <linear unit> {,<axis>, <axis>, <axis>} {,<authority>}]
    
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end

    case 'GEOGCS'       % ["<name>", <datum>, <prime meridian>, <angular unit> {,<twin axes>} {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'INVERSE_MT'
    case 'LOCAL_DATUM'  % ["<name>", <datum type> {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).datumType = A.elements{2};
        for ii = 3:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'LOCAL_CS'     % ["<name>", <local datum>, <unit>, <axis>, {,<axis>}* {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'PARAMETER'
        
        name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(name) = A.elements{2};
        
    case 'PARAM_MT'
    case 'PASSTHROUGH_MT'
    case 'PRIMEM'       % ["<name>", <longitude> {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).longitude = str2num(A.elements{2});
        for ii = 3:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'PROJCS'       % ["<name>", <geographic cs>, <projection>, {<parameter>,}* <linear unit> {,<twin axes>}{,<authority>}]
    
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
    
    case 'PROJECTION'   % ["<name>" {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'SPHEROID'     % ["<name>", <semi-major axis>, <inverse flattening> {,<authority>}]
    
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).semiMajorAxis = str2num(A.elements{2});
        S.(A.key).inverseFlattening = str2num(A.elements{3});
        for ii = 4:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
    
    case 'TOWGS84'      % [<seven param>]
        
        S.(A.key).dx  = str2num(A.elements{1});
        S.(A.key).dy  = str2num(A.elements{2});
        S.(A.key).dz  = str2num(A.elements{3});
        S.(A.key).ex  = str2num(A.elements{4});
        S.(A.key).ey  = str2num(A.elements{5});
        S.(A.key).ez  = str2num(A.elements{6});
        S.(A.key).ppm = str2num(A.elements{7});
        
    case 'UNIT'         % ["<name>", <conversion factor> {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).conversionFactor = str2num(A.elements{2});
        for ii = 3:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
   
    case 'VERT_DATUM'   % ["<name>", <datum type> {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        S.(A.key).datumType = A.elements{2};
        for ii = 3:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
        
    case 'VERT_CS'      % ["<name>", <vert datum>, <linear unit>, {<axis>,} {,<authority>}]
        
        S.(A.key).name = regexprep(A.elements{1}, '(^\")|(\"$)', '');
        for ii = 2:length(A.elements)
            S.(A.key) = addstructure(S.(A.key), A.elements{ii});
        end
    
    otherwise
        error('Unexpected flag');
end
   

% Subfunction: parse a field that is itself a structure

function St = addstructure(St, str)
Tmp = parsewkt(str);
St = mergestruct(St, Tmp);
    
% Subfunction: separate comma-delimited elements within outer set of
% brackets 
        
function A = parsegroup(str)

A = regexp(str, '(?<key>[A-Z_0-9]*)\[(?<val>.*)\]$', 'names');
   
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

if length(A) > 0 && ismember(A.key, keywords)
    openbrack = 0;
    idx = 1;
    A.elements = cell(0);
    for ic = 1:length(A.val)
        switch A.val(ic)
            case '['
                openbrack = openbrack + 1;
            case ']'
                openbrack = openbrack - 1;
            case ','
                if openbrack == 0
                    A.elements = [A.elements; {A.val(idx:ic-1)}];
                    idx = ic;
                end
        end        
    end
    A.elements = [A.elements; {A.val(idx:end)}];

    A.elements = regexprep(A.elements, '^,', '');
else
    blah
end

    







