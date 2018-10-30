cmake_minimum_required (VERSION 2.8)
project (NAD_program)
enable_language (Fortran)

###################################################################################################
# Options.                                                                                        #
  if(NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
         "Choose the type of build, options are: Release, Debug and DebugAll ..."
         FORCE)
  endif()
  set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug DebugAll)


  set(LINALG MKL CACHE STRING "Linear algebra library, at the moment only MKL is supported")
  set_property(CACHE LINALG PROPERTY STRINGS MKL)

  option(OMP "Use OMP." ON)
  option(STATIC "Create a static executable." OFF)

  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
  set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)

  include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})
  link_directories(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})


###################################################################################################
# Compiler specific options:                                                                      #
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    add_compile_options(-fpp)
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG    "-O2 -traceback -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -traceback -g -check all -debug all")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-cpp)
    add_compile_options(-ffree-line-length-0)
    set(CMAKE_Fortran_FLAGS_RELEASE  "-O3")
    set(CMAKE_Fortran_FLAGS_DEBUG "-O2 -fbacktrace -g")
    set(CMAKE_Fortran_FLAGS_DEBUGALL "-O0 -fbacktrace -Wall -pedantic -fcheck=all -fbounds-check -g")
  endif()


###################################################################################################
# Libraries.                                                                                      #
  add_definitions(-DLINALG)
  # CMake FindBLAS/FindLAPACK modules don't work for ilp64 version of MKL or for MKL F95 interface.
  # Set everything manually.
  if(LINALG STREQUAL "MKL")
    set (MKLROOT "" CACHE PATH "MKL root directory.")
    set (F95ROOT "" CACHE PATH "MKL F95 interface root directory.")
    if (MKLROOT STREQUAL "")
      set (MKLROOT $ENV{MKLROOT})
      if (NOT MKLROOT)
        message (FATAL_ERROR 
                 "When using MKL, set MKLROOT environment variable "
                 "or call cmake with -DMKLROOT=/path/to/mklroot_dir")
      endif()
    endif()
    if (F95ROOT STREQUAL "")
      set (F95ROOT $ENV{F95ROOT})
      if (NOT F95ROOT)
        message (FATAL_ERROR 
                 "When using MKL, set F95ROOT environment variable "
                 "or call cmake with -DF95ROOT=/path/to/f95root_dir")
      endif()
    endif()
    message (STATUS "MKLROOT = ${MKLROOT}")
    message (STATUS "F95ROOT = ${F95ROOT}")

    # Compile flags:
    # F95 Interface:
    include_directories(${F95ROOT}/include/intel64/ilp64)
    set(F95_LIBS ${F95ROOT}/lib/intel64/libmkl_blas95_ilp64.a)
    set(F95_LIBS ${F95_LIBS} ${F95ROOT}/lib/intel64/libmkl_lapack95_ilp64.a)
    set(LINKLIBS ${LINKLIBS} ${F95_LIBS})
    # MKL libraries:
    include_directories(${MKLROOT}/include) 
    set(MKL_LIBS "-Wl,--no-as-needed -Wl,--start-group")
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(MKL_LIBS ${MKL_LIBS} "-lmkl_gf_ilp64")
      add_compile_options(-fdefault-integer-8 -m64)
    else()
      add_compile_options(-i8 -m64)
      set(MKL_LIBS ${MKL_LIBS} "-lmkl_intel_ilp64")
    endif()
    set(MKL_LIBS ${MKL_LIBS} "-lmkl_sequential -lmkl_core -Wl,--end-group -lpthread -lm -ldl")
    set(LINKLIBS ${LINKLIBS} ${MKL_LIBS})
  endif()

  if(OMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    message(STATUS "OpenMP libraries:  ${OpenMP_Fortran_LIB_NAMES}")
  endif()

  if(STATIC)
    set(LINKLIBS ${LINKLIBS} "-static")
  endif()


################################################################################
# Build                                                                        #
  add_subdirectory(theochemIO)
  set(LINKLIBS ${LINKLIBS} "-ltheochemio") 

  file(GLOB SOURCES "src/core/*f90" "src/utils/*f90")
  add_executable(cis_nto.exe "src/prog_nto.f90" ${SOURCES})
  add_executable(cis_overlap.exe "src/prog_overlap.f90" ${SOURCES})
  add_executable(cis_dyson.exe "src/prog_dyson.f90" ${SOURCES})

  add_dependencies(cis_nto.exe theochemio)
  add_dependencies(cis_overlap.exe theochemio)
  add_dependencies(cis_dyson.exe theochemio)

  target_link_libraries(cis_nto.exe ${LINKLIBS})
  target_link_libraries(cis_overlap.exe ${LINKLIBS})
  target_link_libraries(cis_dyson.exe ${LINKLIBS})
