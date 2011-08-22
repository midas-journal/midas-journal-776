

# CMakeLists.txt for extended 2D/3D-registration framework
# Author: Philipp Steininger
# 2009-12-17


SET(OWN_LIB_NAME
  2D3DRegistrationFrameworkExtended
)

SET(OWN_LIB_SRCS
  itkMattesMutualInformationImageToImageMetricComplete.h
  itkMattesMutualInformationImageToImageMetricComplete.txx
  itkMultiResolutionImage2D3DRegistrationMethodCustom.h
  itkMultiResolutionImage2D3DRegistrationMethodCustom.txx
  itkMultiResolutionRegistrationCommandCustom.h
  itkMultiResolutionRegistrationCommandCustom.txx
)

# convince MSVC that txx files are also C++ files!
IF(MSVC)
  SET_SOURCE_FILES_PROPERTIES(
    itkMattesMutualInformationImageToImageMetricComplete.txx 
    itkMultiResolutionImage2D3DRegistrationMethodCustom.txx 
    itkMultiResolutionRegistrationCommandCustom.txx 
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

TARGET_LINK_LIBRARIES(${OWN_LIB_NAME} ITKCommon ITKNumerics ITKBasicFilters ITKAlgorithms DRR 2D3DRegistrationFramework)