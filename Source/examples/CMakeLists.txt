
# Projection generator (copy of InsightApplications)
ADD_EXECUTABLE(GenerateProjection
 
  GenerateProjection.cxx
)

TARGET_LINK_LIBRARIES(GenerateProjection 
  ITKCommon
  ITKBasicFilters
  ITKIO
)


# Intensity-based 2D/3D-registration (modified copy of InsightApplications)
INCLUDE_DIRECTORIES(
  ${SRCMetric_SOURCE_DIR}
)
ADD_EXECUTABLE(IntensityBased2D3DRegistration2
 
  IntensityBased2D3DRegistration2.cxx
)

TARGET_LINK_LIBRARIES(IntensityBased2D3DRegistration2 
  ITKCommon
  ITKBasicFilters
  ITKNumerics
  ITKStatistics
  ITKIO
  Metrics
)


# Extended ITK-based 2D/3D-registration framework, see
# http://ibia.umit.at/ResearchGroup/Phil/web/Simple2D3DRegistrationFramework.html
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
)
ADD_SUBDIRECTORY(fwex)
LINK_DIRECTORIES(fwex)


# 2D/3D-registration application based on extended ITK-based framework (above)
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${CMAKE_CURRENT_SOURCE_DIR}/fwex
${CMAKE_CURRENT_SOURCE_DIR}/fwex/2D3DRegistrationFramework
${CMAKE_CURRENT_SOURCE_DIR}/fwex/2D3DRegistrationFrameworkExtended
${CMAKE_CURRENT_SOURCE_DIR}/fwex/DRR
) 

SET(2D3DRegistration2_SRCS
  2D3DRegistration2.cxx
  oraMultiResolutionImage2D3DRegistrationMethodWithSRC.h
  oraMultiResolutionImage2D3DRegistrationMethodWithSRC.txx
  oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.h
  oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.txx
  oraMultiResolutionRegistrationCommandWithSRC.h
  oraMultiResolutionRegistrationCommandWithSRC.txx
  oraOptimizerIterationCommandWithSRC.h
  oraOptimizerIterationCommandWithSRC.txx
)

# convince MSVC that txx files are also C++ files!
IF(MSVC)
  SET_SOURCE_FILES_PROPERTIES(
    oraMultiResolutionImage2D3DRegistrationMethodWithSRC.txx
    oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.txx
    oraMultiResolutionRegistrationCommandWithSRC.txx
    oraOptimizerIterationCommandWithSRC.txx
    PROPERTIES 
    LANGUAGE "CXX" 
    HEADER_FILE_ONLY FALSE
    COMPILE_FLAGS ""
  )
ENDIF(MSVC)

ADD_EXECUTABLE(2D3DRegistration2
	${2D3DRegistration2_SRCS}
)
TARGET_LINK_LIBRARIES(2D3DRegistration2
	ITKCommon 
	ITKNumerics 
	ITKIO 
	ITKStatistics
	2D3DRegistrationFramework 
	2D3DRegistrationFrameworkExtended
)
