enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)
cmake_minimum_required (VERSION 3.18)
project (NAD_program)
enable_language (Fortran)
enable_language (C)

###################################################################################################
# Options.                                                                                        #
  if(NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Release, Debug and DebugAll ..."
        FORCE)
  endif()
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug DebugAll)

  set(LINALG LAPACK CACHE STRING "Linear algebra library.")
  set_property(CACHE LINALG PROPERTY STRINGS LAPACK MKL)

  option(OMP "Use OMP." ON)
  option(STATIC "Create a static executable." OFF)

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
include_directories(${CMAKE_SOURCE_DIR}/src/easywfo/eigen)
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


###################################################################################################
# Compiler specific options:                                                                      #
  if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-fpp")
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG    "-O2 -traceback -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -traceback -g -check all -debug all")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(dialect "-cpp -ffree-line-length-0")
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -fbacktrace -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -fbacktrace -Wall -pedantic -fcheck=all -fbounds-check -g")
  endif()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


# NOT TESTED FOR INTEL!!!!
  if(CMAKE_C_COMPILER_ID MATCHES "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  elseif(CMAKE_C_COMPILER_ID STREQUAL "GCC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  endif()



###################################################################################################
# Libraries.                                                                                      #
  add_definitions(-DLINALG=1)
  if(LINALG STREQUAL "LAPACK")
    if (BLAS_64)
      set(BLA_SIZEOF_INTEGER 8)
      add_definitions(-DBLA_INT64=1)
    else()
      set(BLA_SIZEOF_INTEGER 4)
    endif()
    find_package(LAPACK REQUIRED)
    add_definitions(-DOPENBLAS=1)
    find_path(CBLAS_INCLUDE_DIR cblas.h PATHS /usr/include /usr/include/x86_64-linux-gnu /usr/include/openblas)
    if(CBLAS_INCLUDE_DIR)
        include_directories(${CBLAS_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src/easywfo/eigen)
    endif()
    find_library(CBLAS_LIBRARY NAMES cblas openblas)
    set(LINALG_LIBS LAPACK::LAPACK ${CBLAS_LIBRARY})
  elseif(LINALG STREQUAL "MKL")
    if (BLAS_64)
      set(MKL_INTERFACE "ilp64")
      add_definitions(-DBLA_INT64=1)
    else()
      set(MKL_INTERFACE "lp64")
    endif()
    find_package(MKL CONFIG REQUIRED)
    set(LINALG_LIBS MKL::MKL)
    add_definitions(-DMKL=1)
  endif()


  if(OMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "OpenMP libraries:  ${OpenMP_Fortran_LIB_NAMES}")
  endif()

  if(STATIC)
    set(LINKLIBS ${LINKLIBS} "-static")
  endif()


################################################################################
# Build                                                                        #

  file(GLOB SOURCES_UTIL "src/util_chem/*/*f90" "src/util_fort/*f90")
  add_library(theochemutil STATIC ${SOURCES_UTIL})
  target_link_libraries(theochemutil ${LINALG_LIBS})

  file(GLOB SOURCES_OLAP
    "src/overlap/*f90"
    "src/util_chem/l2m/*c"
    "src/util_chem/l2m/*f90"
    "src/util_chem/cis/cis_overlap_l2m.f90"
)
  file(GLOB SOURCES_CONV "src/convert/*f90")
  file(GLOB SOURCES_GNTO "src/nto/*f90")

  add_executable(cis_overlap.exe src/easywfo/easywfo_bridge.cpp src/easywfo/easywfo_mod.f90 ${SOURCES_OLAP})
  target_compile_definitions(cis_overlap.exe PRIVATE L2M=1)
target_compile_options(cis_overlap.exe PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-cpp>)
  target_link_libraries(cis_overlap.exe stdc++ ${LINALG_LIBS} theochemutil ${LINKLIBS})

  add_executable(cis_nto.exe ${SOURCES_GNTO})
  target_link_libraries(cis_nto.exe ${LINALG_LIBS} theochemutil ${LINKLIBS})

  add_executable(convert.exe ${SOURCES_CONV})
  target_link_libraries(convert.exe ${LINALG_LIBS} theochemutil ${LINKLIBS})

################################################################################
# Check                                                                        #
  file(GLOB SOURCES_TEST "test/unit_test/*f90")
  add_custom_target(check COMMAND
                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cis_test.exe)
  add_executable(cis_test.exe EXCLUDE_FROM_ALL ${SOURCES_TEST} ${SOURCES_OVLC})
  add_dependencies(check cis_test.exe)
  target_link_libraries(cis_test.exe theochemutil ${LINKLIBS})
