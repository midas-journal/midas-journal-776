
SET(Metrics_SRCS
  oraAverageRanksImageToImageFilter.h
  oraAverageRanksImageToImageFilter.txx
  oraStochasticRankCorrelationImageToImageMetric.h
  oraStochasticRankCorrelationImageToImageMetric.txx
)

# convince MSVC that txx files are also C++ files!
IF(MSVC)
  SET_SOURCE_FILES_PROPERTIES(
    oraAverageRanksImageToImageFilter.txx
    oraStochasticRankCorrelationImageToImageMetric.txx
    PROPERTIES 
    LANGUAGE "CXX" 
    HEADER_FILE_ONLY FALSE
    COMPILE_FLAGS ""
  )
ENDIF(MSVC)

ADD_LIBRARY(Metrics
  ${Metrics_SRCS}
)

SET_TARGET_PROPERTIES(Metrics PROPERTIES LINKER_LANGUAGE CXX)


TARGET_LINK_LIBRARIES(Metrics 
  ITKAlgorithms
  ITKBasicFilters
  ITKStatistics
)

# UNIT TEST
IF (BUILD_TESTING)
  
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
  ADD_SUBDIRECTORY(testing)
  LINK_DIRECTORIES(testing)
    
ENDIF (BUILD_TESTING)  


