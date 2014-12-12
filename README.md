## shapeprjread.m Documentation

This function reads data from a shapefile, incorporating information from a projection file (.prj component of the shapefile).  

If the projection file indicates a geographic
  coordinate system, data will be read into a geographic data structure.
  If the projection file indicates a projected coordinate system, the
  original coordinates (X,Y) will be read and then reverse projected to
  Lat/Lon coordinates.

### Syntax

```
Shp = shapeprjread(file)
Shp = shapeprjread(file, p1, v1, ...)
[Shp, m] = shapeprjread(file, p1, v1, ...)
```
See function help for description of input and output variables.  Syntax is identical to shaperead.m.

