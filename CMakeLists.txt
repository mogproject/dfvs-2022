# OPTIL.io uses CMake version 3.5.1, as of April 11, 2021.
# see https://www.optil.io/optilion/environment
cmake_minimum_required(VERSION 3.5)

# options
option(BUILD_TESTS "Build test programs" OFF)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# specify compilers on local machine
if (EXISTS /usr/local/bin/gcc-11)
  set(CMAKE_C_COMPILER /usr/local/bin/gcc-11)
endif ()

if (EXISTS /usr/local/bin/c++-11)
  set(CMAKE_CXX_COMPILER /usr/local/bin/c++-11)
endif ()

# get external source files
file(
  GLOB MIS_SRC
  src/pace-2019/lib/*/*.h
  src/pace-2019/lib/*/*.cpp
  src/pace-2019/lib/*/*/*.h
  src/pace-2019/lib/*/*/*.cpp
  src/pace-2019/extern/KaHIP/lib/tools/random_functions.h
  src/pace-2019/extern/KaHIP/lib/tools/random_functions.cpp
)

include_directories(
  src/pace-2019/extern/KaHIP/interface
  src/pace-2019/extern/KaHIP/lib
  src/pace-2019/extern/KaHIP/lib/data_structure
  src/pace-2019/extern/KaHIP/lib/io
  src/pace-2019/extern/KaHIP/lib/tools
  src/pace-2019/extern/cliqueSolver
  src/pace-2019/lib
  src/pace-2019/lib/tools
  src/pace-2019/lib/mis
  src/pace-2019/lib/mis/initial_mis
  src/pace-2019/lib/mis/evolutionary
  src/pace-2019/lib/mis/ils
  src/pace-2019/extern/KaHIP/lib/partition
)

# get all source files
file(GLOB_RECURSE DFVS_SRC src/main/solver.cpp src/main/*/*.[ch]pp)

# project
project(DFVS)
add_executable(dfvs ${DFVS_SRC} ${MIS_SRC})

set(CMAKE_CXX_FLAGS "-Wall -funroll-loops -fno-stack-limit -O3")
set (CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

if (NOT DEFINED BUILD_DEV)
  set(BUILD_DEV false)
endif ()

if (NOT DEFINED PROFILE_ON)
  set(PROFILE_ON false)
endif ()

if (NOT DEFINED TRACE_ON)
  set(TRACE_ON false)
endif ()

add_compile_definitions(PROFILE=${PROFILE_ON})
add_compile_definitions(TRACE_ON=${TRACE_ON})
add_compile_definitions(BUILD_DEV=${BUILD_DEV})

# tests with GoogleTest
if (BUILD_TESTS)
  add_subdirectory(test)
endif (BUILD_TESTS)
