cmake_minimum_required (VERSION 3.14)

project (bmw C Fortran)

find_package (OpenMP)

add_library (bmw STATIC)

target_link_libraries (bmw PUBLIC stell)

if (OpenMP_Fortran_FOUND)
    target_link_libraries (bmw PUBLIC OpenMP::OpenMP_Fortran)
endif ()

add_executable (xbmw)
target_link_libraries (xbmw PUBLIC bmw)

add_subdirectory (Sources)
