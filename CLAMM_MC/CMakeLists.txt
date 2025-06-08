cmake_minimum_required(VERSION 3.10)
project(CLAMM_MC)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Force a Release build by default (and record it in the cache)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      Release                     # default
      CACHE STRING "Build type" FORCE)
endif()

# Define output directories
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)

# Compiler flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
   CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG   "-g")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
  set(CMAKE_CXX_FLAGS_DEBUG   "/Zi")
endif()

# Find Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Gather sources
file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)

# Build executable
add_executable(CLAMM_MC ${SOURCES})
if(Boost_FOUND)
  target_link_libraries(CLAMM_MC ${Boost_LIBRARIES})
endif()
include_directories(${CMAKE_SOURCE_DIR}/include)

# Report
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
