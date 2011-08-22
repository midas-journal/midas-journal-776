

# CMakeLists.txt for DRR capabilities (programmed by Markus Neuner)
# Author: Philipp Steininger
# 2009-12-09


SET(OWN_LIB_NAME
  DRR
)

SET(OWN_LIB_SRCS
  itkPerspectiveProjectionImageFilter.h
  itkPerspectiveProjectionImageFilter.txx
  itkPerspectiveProjectionInterpolateImageFunction.h
  itkPerspectiveProjectionInterpolateImageFunction.txx
  itkRayCastPerspectiveProjectionImageFilter.h
  itkRayCastPerspectiveProjectionImageFilter.txx
  itkSplatPerspectiveProjectionImageFilter.h
  itkSplatPerspectiveProjectionImageFilter.txx
  itkUpdateCommand.h 
)

# convince MSVC that txx files are also C++ files!
IF(MSVC)
  SET_SOURCE_FILES_PROPERTIES(
    itkPerspectiveProjectionImageFilter.txx 
    itkPerspectiveProjectionInterpolateImageFunction.txx 
    itkRayCastPerspectiveProjectionImageFilter.txx 
    itkSplatPerspectiveProjectionImageFilter.txx 
    PROPERTIES 
    LANGUAGE "CXX" 
    HEADER_FILE_ONLY FALSE
    COMPILE_FLAGS ""
  )
ENDIF(MSVC)

ADD_LIBRARY(${OWN_LIB_NAME}
  ${OWN_LIB_SRCS}
)

SET_TARGET_PROPERTIES(${OWN_LIB_NAME} PROPERTIES LINKER_LANGUAGE CXX)

TARGET_LINK_LIBRARIES(${OWN_LIB_NAME} ITKCommon ITKNumerics ITKBasicFilters)