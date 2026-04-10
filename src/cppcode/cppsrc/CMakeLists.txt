cmake_minimum_required(VERSION 3.20)
project(mpicpp NONE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
# -------------------------------------------------------------------
# 1. Output Directories
# -------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../shared_lib")

# Disable Kokkos warning about not supporting C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# -------------------------------------------------------------------
# 2. Find Dependencies & Set Compilers BEFORE defining targets
# -------------------------------------------------------------------
# Get Trilinos
if(DEFINED ENV{Trilinos_DIR} AND NOT "$ENV{Trilinos_DIR}" STREQUAL "")
  set(Trilinos_DIR "$ENV{Trilinos_DIR}")
endif()

find_package(Trilinos REQUIRED COMPONENTS Amesos2 Teuchos Tpetra Belos ShyLU_DDFROSch Xpetra)

# Find MPI
find_package(MPI REQUIRED)

# Make sure to use same compilers and flags as Trilinos
set(CMAKE_CXX_COMPILER ${Trilinos_CXX_COMPILER})
set(CMAKE_C_COMPILER ${Trilinos_C_COMPILER})
set(CMAKE_Fortran_COMPILER ${Trilinos_Fortran_COMPILER})
set(CMAKE_CXX_FLAGS  "${Trilinos_CXX_COMPILER_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_FLAGS  "${Trilinos_C_COMPILER_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_Fortran_FLAGS  "${Trilinos_Fortran_COMPILER_FLAGS} ${CMAKE_Fortran_FLAGS}")

# Now enable the languages using the Trilinos compilers
enable_language(C)
enable_language(CXX)
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran)
endif()

# Set the C++ standard and require it
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Echo trilinos build info
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
MESSAGE("End of Trilinos details\n")

# -------------------------------------------------------------------
# 3. Dynamic Path Extraction
# -------------------------------------------------------------------
# Resolve the Julia project root relative to this CMake file instead of
# relying on DrWatson activation from within package-internal source.
get_filename_component(PROJECT_DIR "${CMAKE_SOURCE_DIR}/../../.." ABSOLUTE)

execute_process(
  COMMAND julia --project=${PROJECT_DIR}/ -e "using libcxxwrap_julia_jll; print(libcxxwrap_julia_jll.artifact_dir)"
  OUTPUT_VARIABLE LIBCXXWRAP_ARTIFACT_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND julia --project=${PROJECT_DIR}/ -e "print(normpath(Sys.BINDIR, \"..\", \"include\", \"julia\"))"
  OUTPUT_VARIABLE JULIA_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
  COMMAND julia --project=${PROJECT_DIR}/ -e "using Libdl; print(Libdl.dlpath(\"libjulia\"))"
  OUTPUT_VARIABLE JULIA_LIBJULIA
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

get_filename_component(JULIA_LIBRARY_DIR "${JULIA_LIBJULIA}" DIRECTORY)

# Safely find the CxxWrap library regardless of OS (.so, .dylib, .dll)
find_library(CXXWRAP_JULIA_LIB 
  NAMES cxxwrap_julia 
  PATHS "${LIBCXXWRAP_ARTIFACT_DIR}/lib" 
  NO_DEFAULT_PATH
)

# -------------------------------------------------------------------
# 4. Target Definition and Linking
# -------------------------------------------------------------------
# Create the target NOW that languages are enabled
add_library(mpicpp SHARED mpicpp.cpp)

# Linking all directories
target_include_directories(mpicpp PRIVATE
  ${LIBCXXWRAP_ARTIFACT_DIR}/include
  ${JULIA_INCLUDE_DIR}
  ${MPI_CXX_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${Trilinos_INCLUDE_DIRS}
  ${Trilinos_TPL_INCLUDE_DIRS}
)

# Link with the necessary libraries
target_link_libraries(mpicpp PRIVATE 
  ${CXXWRAP_JULIA_LIB} 
  ${JULIA_LIBJULIA}
  ${Trilinos_LIBRARIES} 
  ${Trilinos_TPL_LIBRARIES}
)

set_target_properties(mpicpp PROPERTIES
  BUILD_RPATH "${JULIA_LIBRARY_DIR}"
  INSTALL_RPATH "${JULIA_LIBRARY_DIR}"
)

# MPI Related compilation options
target_compile_definitions(mpicpp PRIVATE ${MPI_CXX_COMPILE_DEFINITIONS})
target_compile_options(mpicpp PRIVATE ${MPI_CXX_COMPILE_OPTIONS})
target_link_options(mpicpp PRIVATE ${MPI_CXX_LINK_FLAGS})