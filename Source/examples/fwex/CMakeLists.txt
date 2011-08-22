
# CMakeLists.txt for 2D3DRegistration
# Author: Philipp Steininger
# 2009-12-09


SET(PROJ_NAME 2D3DRegistration)

PROJECT(${PROJ_NAME})


# Need the Insight Segmentation and Registration Toolkit (ITK) to compile.
FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
  "Cannot build ${PROJ_NAME} without ITK. Please set ITK_DIR.")
ENDIF(ITK_FOUND)


# sub-grouped libraries (compile them separately into libraries):

# (A) DRR library (based on the implementation of Markus Neuner):
SET(DRR_SUBDIR ${CMAKE_CURRENT_SOURCE_DIR}/DRR)
ADD_SUBDIRECTORY(${DRR_SUBDIR})
INCLUDE_DIRECTORIES(${DRR_SUBDIR})
#link_directories (${DRR_SUBDIR})

# (B) 2D3DRegistrationFramework:
SET(FW_SUBDIR ${CMAKE_CURRENT_SOURCE_DIR}/2D3DRegistrationFramework)
ADD_SUBDIRECTORY(${FW_SUBDIR})
INCLUDE_DIRECTORIES(${FW_SUBDIR})
#link_directories (${FW_SUBDIR})

# (C) 2D3DRegistrationFrameworkExtended:
SET(FWC_SUBDIR ${CMAKE_CURRENT_SOURCE_DIR}/2D3DRegistrationFrameworkExtended)
ADD_SUBDIRECTORY(${FWC_SUBDIR})
INCLUDE_DIRECTORIES(${FWC_SUBDIR})
#link_directories (${FWC_SUBDIR})


# libraries that need to be linked:
# - ITK
SET(ITK_SEL_LIBS ITKCommon ITKNumerics ITKIO ITKStatistics)
# - Registration Framework, Registration Framework Extended, DRR
SET(OWN_LIBS DRR 2D3DRegistrationFramework 2D3DRegistrationFrameworkExtended)


# create the binary application:
SET(EXT_LIBS ${ITK_SEL_LIBS} ${OWN_LIBS})
SET(BINARY_NAME ${PROJ_NAME})
SET(BINARY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/2D3DRegistration.cxx)
ADD_EXECUTABLE(${BINARY_NAME} ${BINARY_SOURCES})
TARGET_LINK_LIBRARIES(${BINARY_NAME} ${EXT_LIBS})

# FINISHED

