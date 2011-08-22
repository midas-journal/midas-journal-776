

# CMakeLists.txt for basic 2D/3D-registration framework
# Author: Philipp Steininger
# 2009-12-09


SET(OWN_LIB_NAME
  2D3DRegistrationFramework
)

SET(OWN_LIB_SRCS
  itkMultiResolutionImage2D3DRegistrationMethod.h  
   itkMultiResolutionImage2D3DRegistrationMethod.txx
   itkMultiResolutionRegistrationCommand.h
   itkMultiResolutionRegistrationCommand.txx
   itkOptimizerIterationCommand.h
   itkOptimizerIterationCommand.txx
   itkMultiResolutionImage2D3DRegistrationMethodCommand.h
   itkMultiResolutionImage2D3DRegistrationMethodCommand.txx
   itkExhaustiveOptimizerComplete.hxx
)

# convince MSVC that txx and hxx files are also C++ files!
IF(MSVC)
  SET_SOURCE_FILES_PROPERTIES(
    itkMultiResolutionImage2D3DRegistrationMethod.txx 
    itkMultiResolutionRegistrationCommand.txx 
    itkOptimizerIterationCommand.txx 
    itkMultiResolutionImage2D3DRegistrationMethodCommand.txx 
    itkExhaustiveOptimizerComplete.hxx 
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

TARGET_LINK_LIBRARIES(${OWN_LIB_NAME} ITKCommon ITKNumerics ITKBasicFilters ITKAlgorithms DRR)
