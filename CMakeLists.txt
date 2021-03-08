cmake_minimum_required (VERSION 3.14)

project (bmw C Fortran)

find_package (OpenMP)

add_library (bmw STATIC)
add_dependencies (bmw stell)

target_link_libraries (bmw stell)
target_include_directories (bmw PUBLIC $<TARGET_PROPERTY:stell,BINARY_DIR>)

if (OpenMP_Fortran_FOUND)
    target_link_libraries (bmw OpenMP::OpenMP_Fortran)
endif ()

add_executable (xbmw)
add_dependencies (xbmw bmw)
target_link_libraries (xbmw bmw)

add_subdirectory (Sources)
