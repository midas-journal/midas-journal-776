
cmake_minimum_required(VERSION 2.4)

IF(COMMAND cmake_policy)	
  IF("${CMAKE_VERSION}" VERSION_GREATER "2.6" OR "${CMAKE_VERSION}" VERSION_EQUAL "2.6")
  	CMAKE_POLICY(SET CMP0003 NEW)
  ENDIF("${CMAKE_VERSION}" VERSION_GREATER "2.6" OR "${CMAKE_VERSION}" VERSION_EQUAL "2.6")
  IF("${CMAKE_VERSION}" VERSION_GREATER "2.8") 
  	CMAKE_POLICY(SET CMP0015 NEW) # NOT AVAILABLE IN 2.8.0 !!!
  ENDIF("${CMAKE_VERSION}" VERSION_GREATER "2.8")
ENDIF(COMMAND cmake_policy)

PROJECT(SRCMetric)

#include some macros from another file...
INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

OPTION(BUILD_TESTING "Build tests for software-verification." ON)
OPTION(BUILD_EXAMPLES "Build examples for demonstration purposes." ON)

SET(Required_Packages
  ITK
)

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
  ITKCommon
  ITKBasicFilters
  ITKAlgorithms
  ITKNumerics
  ITKIO
)

# A nice option automatically adding current source/binary dir to each CMake-list
SET(CMAKE_INCLUDE_CURRENT_DIR ON)

SET(SRCMetric_SOURCE_DIR ${CMAKE_SOURCE_DIR}/metric)
# Metric library and more
ADD_SUBDIRECTORY(${SRCMetric_SOURCE_DIR})
LINK_DIRECTORIES(${SRCMetric_SOURCE_DIR})

# EXAMPLES
IF (BUILD_EXAMPLES)
  
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${SRCMetric_SOURCE_DIR})
  ADD_SUBDIRECTORY(examples)
  LINK_DIRECTORIES(examples)
  
ENDIF (BUILD_EXAMPLES)
