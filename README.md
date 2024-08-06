
# shapeprjread.m: Read a projected shapefile


Author: Kelly Kearney


This repository includes the code for the `shapeprjread.m` Matlab function, along with all dependent functions required to run it.


The `shaperead` function included with Matlab's Mapping Toolbox can read most shapefiles into a geographic data structure (i.e. mapstruct or geostruct). However, that function ignores any projection information included in the shapefile (i.e., the .prj file). The `shapeprjread` function provides a wrapper that adds this projection data to the input via the companion `prjread` function and calculates lat/lon coordinates appropriately. The resulting output will be a combination mapstruct/geostruct holding both the projected X/Y coordinates from the original file and the reverse-projected lat/lon coordinates corresponding to those.



## Contents

            
- Getting started        
- Syntax        
- Description        
- Contributions

## Getting started


**Prerequisites**


This function requires the Mapping Toolbox, and should run in any version of Matlab.


**Downloading and installation**


This code can be downloaded from [Github](https://github.com/kakearney/shapeprjread-pkg/)


**Matlab Search Path**


The following folders need to be added to your Matlab Search path (via `addpath`, `pathtool`, etc.):



```matlab
shapeprjread-pkg/shapeprjread
```



## Syntax



```
Shp = shapeprjread(file)
Shp = shapeprjread(file, Name, Value, ...)
[Shp, m, fac] = shapeprjread(file, p1, v1, ...)
```



## Description


`S = shapeprjread(file)` reads the shapefile specified by `filename` and returns a geographic data structure array `S`, which includes both the projected X/Y coordinates from the file and reverse-projected lat/lon coordinates corresponding to those.


`S = shapeprjread(file, Name, Value)` returns a subset of the information in the shapefile. See `shaperead.m` (in Matlab's Mapping Toolbox) for all potential subsetting options. All parameters except 'UseGeoCoords' will be accepted. Note that all selectors, including BoundingBox, will be applied to the data as it is stored in the file, not to the reverse-projected data.


`[S, m, fac] = shapeprjread(...)` returns the map projection structure `m` that is used for the reverse projection and the projection unit conversion factor `fac`. This factor is usually 1 but will vary in files that include a unit conversion factor; this factor can be used to replicate the file's conversion as follows:



  - reverse projection: `[lat,lon] = minvtran(m, x*fac, y*fac);`
  - forward projection: `[x,y] = mfwdtran(lat,lon); x = x/fac; y = y/fac;`


## Contributions


Community contributions to this package are welcome!


To report bugs, please submit [an issue](https://github.com/kakearney/shapeprjread-pkg/issues) on GitHub and include:



  - your operating system
  - your version of Matlab and all relevant toolboxes (type `ver` at the Matlab command line to get this info)
  - code/data to reproduce the error or buggy behavior, and the full text of any error messages received

Please also feel free to submit enhancement requests, or to send pull requests (via GitHub) for bug fixes or new features.


I do monitor the MatlabCentral FileExchange entry for any issues raised in the comments, but would prefer to track issues on GitHub.



<sub>[Published with MATLAB R2024a]("http://www.mathworks.com/products/matlab/")</sub>
