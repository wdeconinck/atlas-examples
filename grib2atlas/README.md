Toy example to convert a GRIB file into a custom Atlas-IO file format.

Why?
----
Atlas-IO format does not need to adhere to any standard as GRIB,
and can encode any kind of data. This could be useful for experimenting
with new grids.

DISCLAIMERS:
------------
- No support is to be expected for this example
- The grids that are encoded are expected to be global Gaussian grids only.

Compilation
-----------

    export atlas_ROOT=<path-to-atlas>
    export eccodes_ROOT=<path-to-eccodes>
    mkdir build; cd build
    cmake <path-to-source>
    make

There should now be an executable "grib2atlas" in the build directory


