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
  pace-2019/lib/*/*.h
  pace-2019/lib/*/*.cpp
  pace-2019/lib/*/*/*.h
  pace-2019/lib/*/*/*.cpp
  pace-2019/extern/KaHIP/lib/tools/random_functions.h
  pace-2019/extern/KaHIP/lib/tools/random_functions.cpp
)

include_directories(
  pace-2019/extern/KaHIP/interface
  pace-2019/extern/KaHIP/lib
  pace-2019/extern/KaHIP/lib/data_structure
  pace-2019/extern/KaHIP/lib/io
  pace-2019/extern/KaHIP/lib/tools
  pace-2019/extern/cliqueSolver
  pace-2019/lib
  pace-2019/lib/tools
  pace-2019/lib/mis
  pace-2019/lib/mis/initial_mis
  pace-2019/lib/mis/evolutionary
  pace-2019/lib/mis/ils
  pace-2019/extern/KaHIP/lib/partition
)

# get all source files
if (BUILD_DEV)
  file(GLOB_RECURSE DFVS_SRC main/solver_dev.cpp main/*/*.[ch]pp)
else ()
  file(GLOB_RECURSE DFVS_SRC main/solver.cpp main/*/*.[ch]pp)
endif ()

file(GLOB_RECURSE VERIFY_SRC
  main/verifier.cpp
  main/algorithm/component.[ch]pp
  main/data/*.[ch]pp
  main/readwrite/*.[ch]pp
  main/util/*.[ch]pp
)

# project
project(DFVS)
add_executable(dfvs ${DFVS_SRC} ${MIS_SRC})

set(CMAKE_CXX_FLAGS "-Wall -funroll-loops -fno-stack-limit -O3")
# set(CMAKE_EXE_LINKER_FLAGS "")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
    set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
else ()
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

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
