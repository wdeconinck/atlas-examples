cmake_minimum_required( VERSION 3.18 )

project( grib2atlas VERSION 0.0.0 LANGUAGES C CXX )

find_package( eccodes REQUIRED )
find_package( atlas REQUIRED )

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable( grib2atlas )
target_sources( grib2atlas PUBLIC main.cc )
target_link_libraries( grib2atlas PUBLIC atlas eccodes )
