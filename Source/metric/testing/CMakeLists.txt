
# average ranks filter test
INCLUDE_DIRECTORIES(
  ${SRCMetric_SOURCE_DIR}
)
ADD_EXECUTABLE(TestAverageRanksFilter
  
  TestAverageRanksFilter.cxx
)
TARGET_LINK_LIBRARIES(TestAverageRanksFilter
  Metrics
  ITKIO
)
ADD_TEST(TestAverageRanksFilter ${EXECUTABLE_OUTPUT_PATH}/TestAverageRanksFilter -v)


# stochastic rank correlation metric test
INCLUDE_DIRECTORIES(
  ${SRCMetric_SOURCE_DIR}
)
ADD_EXECUTABLE(TestStochasticRankCorrelationMetric
  
  TestStochasticRankCorrelationMetric.cxx
)
TARGET_LINK_LIBRARIES(TestStochasticRankCorrelationMetric
  Metrics
  ITKIO
)
ADD_TEST(TestStochasticRankCorrelationMetric ${EXECUTABLE_OUTPUT_PATH}/TestStochasticRankCorrelationMetric -v)

