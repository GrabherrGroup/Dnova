cmake_minimum_required(VERSION 3.3)

# Use the default gcc/g++ version
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

# require at least gcc 5.0
if (CMAKE_COMPILER_IS_GNUCC AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 5.0)
    message(WARNING "Require at least gcc-5.0")
endif()

project(InSeqt)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -O3 -w")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_HOME_DIRECTORY}/bin)

# include directory in find path where all dependency modules exist
include_directories(./)

# DMap binaries
set(SOURCE_FILES_SITELAPS ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/RestSiteAlignUnit.cc src/SiteLaps.cc)  
set(SOURCE_FILES_DMAP ryggrad/src/base/ErrorHandling.cc ryggrad/src/base/FileParser.cc  ryggrad/src/base/StringUtil.cc ryggrad/src/general/DNAVector.cc ryggrad/src/util/mutil.cc src/RestSiteAlignUnit.cc src/DMap.cc)  

add_executable(SiteLaps             ${SOURCE_FILES_SITELAPS}) 
add_executable(DMap                 ${SOURCE_FILES_DMAP}) 

SET_TARGET_PROPERTIES(SiteLaps PROPERTIES COMPILE_FLAGS "-fopenmp" LINK_FLAGS "-fopenmp")
SET_TARGET_PROPERTIES(DMap PROPERTIES COMPILE_FLAGS "-fopenmp" LINK_FLAGS "-fopenmp")
  
