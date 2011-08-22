

#ifndef ORAOPTIMIZERITERATIONCOMMANDWITHSRC_TXX_
#define ORAOPTIMIZERITERATIONCOMMANDWITHSRC_TXX_

#include "oraOptimizerIterationCommandWithSRC.h"

// support for concrete casting:
#include <itkRegularStepGradientDescentOptimizer.h>
#include "oraStochasticRankCorrelationImageToImageMetric.h"

namespace ora /** open radART **/
{


template <typename TFramework>
OptimizerIterationCommandWithSRC<TFramework>
::OptimizerIterationCommandWithSRC()
 : Superclass()
{
  ;
}

template <typename TFramework>
OptimizerIterationCommandWithSRC<TFramework>
::~OptimizerIterationCommandWithSRC()
{ 
  ;
}

template <typename TFramework>
bool 
OptimizerIterationCommandWithSRC<TFramework>
::GetCurrentIterationParameters(BaseOptimizerConstPointer optimizer,
  int &currIt, int &maxIt,
  typename BaseMetricType::MeasureType &thisIterationCostValue,
  typename BaseOptimizerType::ParametersType &pos)
{
  bool validResults = false;

  if (!optimizer)
    return false;

  /* typedefs for concrete casting */
  typedef itk::RegularStepGradientDescentOptimizer RSGDOptimizerType;
  typedef RSGDOptimizerType::ConstPointer RSGDOptimizerConstPointer;

  // do concrete casts (first these are pure dummies, then casted):
  RSGDOptimizerConstPointer rsgdOptimizer =
      itk::SmartPointer<const RSGDOptimizerType>(
      dynamic_cast<const RSGDOptimizerType *>(
      RSGDOptimizerType::New().GetPointer()));

  if (optimizer->GetNameOfClass() == rsgdOptimizer->GetNameOfClass())
  { 
    rsgdOptimizer = itk::SmartPointer<const RSGDOptimizerType>(
      dynamic_cast<const RSGDOptimizerType *>(optimizer.GetPointer()));
    currIt = rsgdOptimizer->GetCurrentIteration() + 1;
    maxIt = rsgdOptimizer->GetNumberOfIterations();
    validResults = true;
    thisIterationCostValue = rsgdOptimizer->GetValue();
    pos = optimizer->GetCurrentPosition();
  }
  else
  {
    validResults = Superclass::GetCurrentIterationParameters(optimizer, currIt,
        maxIt, thisIterationCostValue, pos);
  }

  return validResults;
}

template <typename TFramework>
void
OptimizerIterationCommandWithSRC<TFramework>
::StopOptimizer()
{ 
  if (!this->m_Framework)
    return;
  
  /* typedefs for concrete casting */
  typedef itk::RegularStepGradientDescentOptimizer RSGDOptimizerType;
  typedef RSGDOptimizerType::Pointer RSGDOptimizerPointer;
  
  RSGDOptimizerPointer rsgd = RSGDOptimizerType::New();
  BaseOptimizerPointer opt = this->m_Framework->GetOptimizer();
  
  if (!opt)
    return;
  
  if (opt->GetNameOfClass() == rsgd->GetNameOfClass())
  {
    rsgd = itk::SmartPointer<RSGDOptimizerType>(
      dynamic_cast<RSGDOptimizerType *>(opt.GetPointer()));
    
    rsgd->StopOptimization();
  }
  else
  {
    Superclass::StopOptimizer();
  }
}

template <typename TFramework>
void
OptimizerIterationCommandWithSRC<TFramework>
::WritePDFImage(int currIt, bool costValueChanged)
{
  typedef StochasticRankCorrelationImageToImageMetric<
    typename FrameworkType::InternalImageType,
    typename FrameworkType::InternalImageType,
    typename FrameworkType::RankImageType> SRCMetricType;
  typedef typename SRCMetricType::Pointer SRCMetricPointer;

  SRCMetricPointer src = SRCMetricType::New();
  if (this->m_Framework->GetMetric()->GetNameOfClass() == src->GetNameOfClass())
  {
    // metric-PDF-output each n-th iteration:
    if ((this->m_PDFOutput && (this->m_PDFModulo > 0) &&
         (currIt > 1) && this->m_Framework &&
         (this->m_PDFBaseFileName.length() > 0) &&
         ((currIt % this->m_PDFModulo) == 0)) ||
        (this->m_PDFOutput && this->m_PDFAutoOutput &&
         (currIt > 0) && this->m_Framework &&
         (this->m_PDFBaseFileName.length() > 0) &&
         costValueChanged))
    {
      if (this->m_Verbose)
        std::cout << std::endl << "Writing current sample distribution " <<
          "to disc (iteration=" << currIt << ")." << std::endl << std::endl;

      FrameworkPointer fw = this->m_Framework;

      fw->GetCurrentMetricHistogram(false);

      // write out a comma-separated file:
      src = itk::SmartPointer<SRCMetricType>(
        dynamic_cast<SRCMetricType *>(fw->GetMetric().GetPointer()));
      std::vector<SampleDistributionEntry> *sd = src->GetSampleDistribution();
      std::ofstream csv;
      // -- file name --
      char imgfn[2048];
      std::string s = this->m_PDFBaseFileName;
      std::ostringstream ss("");
      ss << currIt; // insert iteration-number right before last "."
      std::string::size_type last = s.rfind('.');
      if (last != std::string::npos && last > 1)
        s.insert(last, ss.str());
      sprintf(imgfn, s.c_str(),
        (fw->GetRegistration()->GetCurrentLevel() + 1));
      // -- / file name --
      csv.open(imgfn, std::ios::out);
      if (csv.is_open())
      {
        csv << "fixed;moving;fixed-rank;moving-rank\n";
        for (unsigned int i = 0; i < sd->size(); i++)
          csv << (*sd)[i].FixedIntensity << ";" << (*sd)[i].MovingIntensity <<
            ";" << (*sd)[i].FixedRankIntensity << ";" <<
            (*sd)[i].MovingRankIntensity << "\n";
        csv.close();
      }
    }
  }
  else
  {
    Superclass::WritePDFImage(currIt, costValueChanged); // delegate
  }
}


}


#endif /* ORAOPTIMIZERITERATIONCOMMANDWITHSRC_TXX_ */
