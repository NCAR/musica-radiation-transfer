cmake_minimum_required(VERSION 3.17)
project(MUSICA)
enable_language(Fortran)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################
# compiler debug flags

if (DEFINED CMAKE_Fortran_FLAGS_DEBUG_INIT AND
    "${CMAKE_Fortran_FLAGS_DEBUG_INIT}" MATCHES "${CMAKE_Fortran_FLAGS_DEBUG}")
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-g -ggdb -ffpe-trap='invalid,zero,overflow' -finit-real=snan -fcheck=bounds")
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    SET(CMAKE_Fortran_FLAGS_DEBUG "-g -Ktrap=fp -Mbounds -Kieee -traceback")
  endif()
endif()

message("Fortran compiler settings:")
message(STATUS "Fortran Compiler ID: " ${CMAKE_Fortran_COMPILER_ID})
message(STATUS "Fortran Compiler RELEASE Flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
message(STATUS "Fortran Compiler DEBUG Flags: " ${CMAKE_Fortran_FLAGS_DEBUG})

################################################################################
# NetCDF library

find_path(NETCDF_INCLUDE_DIR netcdf.mod NETCDF.mod
  DOC "NetCDF include directory (must contain netcdf.mod"
  PATHS
    $ENV{NETCDF_HOME}/include
    /usr/lib/gfortran/modules
    /usr/lib64/gfortran/modules
    /opt/local/include)
find_library(NETCDF_C_LIB netcdf
  DOC "NetCDF C library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
find_library(NETCDF_FORTRAN_LIB netcdff
  DOC "NetCDF Fortran library"
  PATHS
    $ENV{NETCDF_HOME}/lib
    $ENV{NETCDF_HOME}/lib64
    /opt/local/lib)
  set(NETCDF_LIBS ${NETCDF_C_LIB} ${NETCDF_FORTRAN_LIB})
include_directories(${NETCDF_INCLUDE_DIR})

message(STATUS "Fortran netcdf lib: " ${NETCDF_FORTRAN_LIB})

################################################################################
# copy test input data to build directory

add_custom_target(copy_input_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/photo-utils/test/data ${CMAKE_BINARY_DIR}/test/data)

################################################################################
# musica-photo-utils library

add_subdirectory(photo-utils/src)
include_directories(${CMAKE_BINARY_DIR}/photo-utils/src)

# musica-rad-xfer library

add_subdirectory(src)
include_directories(${CMAKE_BINARY_DIR}/src)

################################################################################
# tests

enable_testing()
add_subdirectory(test)

################################################################################
