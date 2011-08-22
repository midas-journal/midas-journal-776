

#ifndef ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_TXX_
#define ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_TXX_


#include "oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.h"

// support for concrete casting:
#include <itkRegularStepGradientDescentBaseOptimizer.h>
#include "oraStochasticRankCorrelationImageToImageMetric.h"


namespace ora /** open radART **/
{


template <typename TFramework>
MultiResolutionImage2D3DRegistrationMethodCommandWithSRC<TFramework>
::MultiResolutionImage2D3DRegistrationMethodCommandWithSRC()
 : Superclass()
{
  ;
}

template <typename TFramework>
bool 
MultiResolutionImage2D3DRegistrationMethodCommandWithSRC<TFramework>
::GetCurrentResults(typename FrameworkType::Pointer fw,
  typename FrameworkType::BaseOptimizerType::ParametersType &finalPars,
  double &bestValue, unsigned int &numberOfIterations)
{
  bool validResults = false;

  typedef itk::RegularStepGradientDescentOptimizer RSGDOptimizerType;
  typedef RSGDOptimizerType::ConstPointer RSGDOptimizerConstPointer;

  // do concrete casts (first these are pure dummies, then casted):
  RSGDOptimizerConstPointer rsgdOptimizer =
      itk::SmartPointer<const RSGDOptimizerType>(
      dynamic_cast<const RSGDOptimizerType *>(
      RSGDOptimizerType::New().GetPointer()));

  if (fw->GetOptimizer()->GetNameOfClass() == rsgdOptimizer->GetNameOfClass())
  { 
    rsgdOptimizer =
      itk::SmartPointer<const RSGDOptimizerType>(
      dynamic_cast<const RSGDOptimizerType *>(
      fw->GetOptimizer().GetPointer()));

    numberOfIterations = rsgdOptimizer->GetCurrentIteration() + 1;
    finalPars = rsgdOptimizer->GetCurrentPosition();
    bestValue = rsgdOptimizer->GetValue();

    validResults = true;
  }
  else
  {
    validResults = Superclass::GetCurrentResults(fw, finalPars, bestValue,
        numberOfIterations);
  }

  return validResults;
}

template <typename TFramework>
void
MultiResolutionImage2D3DRegistrationMethodCommandWithSRC<TFramework>
::WriteFinalPDFImage(FrameworkPointer fw,
  typename FrameworkType::BaseOptimizerType::ParametersType finalPars)
{
  typedef StochasticRankCorrelationImageToImageMetric<
    typename FrameworkType::InternalImageType,
    typename FrameworkType::InternalImageType,
    typename FrameworkType::RankImageType> SRCMetricType;
  typedef typename SRCMetricType::Pointer SRCMetricPointer;

  SRCMetricPointer src = SRCMetricType::New();
  if (fw->GetMetric()->GetNameOfClass() == src->GetNameOfClass())
  {
    if (this->m_ImageOutput && (this->m_ImageFileNames.size() >= 2) &&
        (this->m_ImageFileNames[1].length() > 0))
    {
      if (this->m_Verbose)
        std::cout << std::endl << "Writing final sample distribution to disc." <<
          std::endl << std::endl;

      // to get the final sample distribution for current level, the
      // transformation must be set to the result transformation parameters:
      typedef typename FrameworkType::BaseTransformType::ParametersType
        ParametersType;
      ParametersType parsStore = fw->GetTransformation()->GetParameters();
      fw->GetTransformation()->SetParameters(finalPars);

      fw->GetCurrentMetricHistogram(true); // update sample distribution

      // write out a comma-separated file:
      src = itk::SmartPointer<SRCMetricType>(
        dynamic_cast<SRCMetricType *>(fw->GetMetric().GetPointer()));
      std::vector<SampleDistributionEntry> *sd = src->GetSampleDistribution();
      std::ofstream csv;
      char fn[2048];
      sprintf(fn, this->m_ImageFileNames[1].c_str(),
        fw->GetRegistration()->GetNumberOfLevels());
      csv.open(fn, std::ios::out);
      if (csv.is_open())
      {
        csv << "fixed;moving;fixed-rank;moving-rank\n";
        for (unsigned int i = 0; i < sd->size(); i++)
          csv << (*sd)[i].FixedIntensity << ";" << (*sd)[i].MovingIntensity <<
            ";" << (*sd)[i].FixedRankIntensity << ";" <<
            (*sd)[i].MovingRankIntensity << "\n";
        csv.close();
      }

      // restore parameters:
      fw->GetTransformation()->SetParameters(parsStore);
    }
  }
  else
  {
    Superclass::WriteFinalPDFImage(fw, finalPars);
  }
}


}


#endif /* ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_TXX_ */
