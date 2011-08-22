

#ifndef ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_TXX_
#define ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_TXX_

#include "oraMultiResolutionRegistrationCommandWithSRC.h"

// support for concrete casting:
#include <itkRegularStepGradientDescentOptimizer.h>
#include "oraStochasticRankCorrelationImageToImageMetric.h"

#include <itkSpatialObjectToImageFilter.h>

namespace ora /** open radART **/
{


template <typename TFramework>
MultiResolutionRegistrationCommandWithSRC<TFramework>
::MultiResolutionRegistrationCommandWithSRC()
 : Superclass()
{
  ;
}

template <typename TFramework>
MultiResolutionRegistrationCommandWithSRC<TFramework>
::~MultiResolutionRegistrationCommandWithSRC()
{
  ;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommandWithSRC<TFramework>
::Execute(itk::Object *object, const itk::EventObject &event)
{
  Superclass::Execute(object, event);

  // cast the basic components:
  RegistrationPointer registration = this->m_Framework->GetRegistration();
  if (registration->GetCurrentLevel() == 0) // for first level ...
  {
    WriteFixedImageMask();
  }
}

template <typename TFramework>
bool
MultiResolutionRegistrationCommandWithSRC<TFramework>
::GetCurrentResults(RegistrationPointer registration,
  typename BaseOptimizerType::ParametersType &finalPars, double &bestValue,
  unsigned int &numberOfIterations)
{
  bool validResults = false;

  /* typedefs for concrete casting */
  typedef itk::RegularStepGradientDescentOptimizer RGSDOptimizerType;
  typedef RGSDOptimizerType::ConstPointer RGSDOptimizerConstPointer;
  typedef StochasticRankCorrelationImageToImageMetric<InternalImageType,
    InternalImageType, RankImageType> SRCMetricType;
  typedef typename SRCMetricType::ConstPointer SRCMetricConstPointer;

  // concrete casts (first dummies, then references):
  RGSDOptimizerConstPointer rgsdOptimizer =
      itk::SmartPointer<const RGSDOptimizerType>(
      dynamic_cast<const RGSDOptimizerType *>(
      RGSDOptimizerType::New().GetPointer()));

  BaseOptimizerPointer optimizer = itk::SmartPointer<BaseOptimizerType>(
    dynamic_cast<BaseOptimizerType *>(registration->GetOptimizer()));

  // final transformation parameters and final/best metric value:
  if (optimizer->GetNameOfClass() == rgsdOptimizer->GetNameOfClass())
  {
    rgsdOptimizer = itk::SmartPointer<RGSDOptimizerType>(
      dynamic_cast<RGSDOptimizerType *>(optimizer.GetPointer()));

    numberOfIterations = rgsdOptimizer->GetCurrentIteration() + 1;
    finalPars = rgsdOptimizer->GetCurrentPosition();
    bestValue = rgsdOptimizer->GetValue();
    validResults = true;
  }
  else
  {
    validResults = Superclass::GetCurrentResults(registration, finalPars,
        bestValue, numberOfIterations);
  }

  return validResults;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommandWithSRC<TFramework>
::WriteFixedImageMask()
{
  // write out fixed image mask, if we have one defined and a file name pattern:
  if (this->m_ImageFileNames.size() >= 7 &&
      this->m_Framework->GetMetric()->GetFixedImageMask())
  {
    RegistrationPointer registration = this->m_Framework->GetRegistration();
    unsigned int level = registration->GetCurrentLevel();

    if (this->m_Verbose)
      std::cout << std::endl << "Writing fixed image mask " <<
        "for level " << level << std::endl << std::endl;

    typedef typename FrameworkType::BaseMetricType MetricType;
    typedef typename MetricType::FixedImageMaskType SpatialObjectType;
    typedef typename FrameworkType::FixedMaskImageType MaskImageType;
    typedef itk::SpatialObjectToImageFilter<SpatialObjectType, MaskImageType>
      SOTIFilterType;
    typedef typename SOTIFilterType::Pointer SOTIFilterPointer;
    typedef itk::ImageFileWriter<MaskImageType> WriterType;
    typedef typename WriterType::Pointer WriterPointer;

    SOTIFilterPointer soti = SOTIFilterType::New();
    SpatialObjectType *mask = static_cast<SpatialObjectType *>(
        (void *)this->m_Framework->GetMetric()->GetFixedImageMask());
    soti->SetInput(mask);
    WriterPointer w = WriterType::New();
    char fn[2048];
    sprintf(fn, this->m_ImageFileNames[6].c_str(), level + 1);
    w->SetFileName(fn);
    w->SetInput(soti->GetOutput());
    try
    {
      soti->UpdateLargestPossibleRegion();
      w->Update();
    }
    catch (itk::ExceptionObject &e)
    {
      std::cerr << "Error during writing fixed image mask (" <<
          fn << ";" << e << ")." << std::endl;
    }
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommandWithSRC<TFramework>
::ReConfigureRegistrationComponents()
{  
  RegistrationPointer registration = this->m_Framework->GetRegistration();
  unsigned int level = registration->GetCurrentLevel(); 
  if (level >= 1)
  {
    // metric
    if (this->m_MetricSels.size() >= level &&
        this->m_MetricConfigs.size() >= this->m_MetricSels.size())
    {
      if (this->m_MetricSels[level - 1] == "SRC")
      {
        if (this->m_Verbose)
          std::cout << std::endl << "Resetting metric component for level " <<
            (level + 1) << ": " << this->m_MetricSels[level - 1];

        typedef StochasticRankCorrelationImageToImageMetric<
          InternalImageType, InternalImageType, RankImageType> SRCMetricType;
        typename SRCMetricType::Pointer src = SRCMetricType::New();
        src->SetFixedImage(registration->GetFixedImage());
        src->SetMovingImage(registration->GetMovingImage());
        src->SetTransform(registration->GetTransform());
        src->SetInterpolator(registration->GetInterpolator());
        typename InternalImageType::RegionType freg = registration->GetMetric()->
          GetFixedImageRegion();
        src->SetFixedImageRegion(freg);
        if (this->m_ImageOutput && (this->m_ImageFileNames.size() >= 2) &&
            (this->m_ImageFileNames[1].length() > 0)) // "PDF"-output
          src->SetExtractSampleDistribution(true);
        registration->SetMetric(src);
        this->m_Framework->SetMetric(registration->GetMetric());
  
        if (this->m_Verbose)
          std::cout << " --> SUCCESS" << std::endl;
      } // else: expected to be processed by super-class
    }
    // this->m_MetricConfigs: processed by superclass ...

    // optimizer
    if (this->m_OptimizerSels.size() >= level &&
        this->m_OptimizerConfigs.size() >= this->m_OptimizerSels.size())
    {
      if (this->m_Verbose)
        std::cout << std::endl << "Resetting optimizer component for level " <<
          (level + 1) << ": " << this->m_OptimizerSels[level - 1];
      if (this->m_OptimizerSels[level - 1] == "RSGD")
      {
        typedef itk::RegularStepGradientDescentOptimizer RGSDOptimizerType;
        RGSDOptimizerType::Pointer rgsd = RGSDOptimizerType::New();
        rgsd->SetCostFunction(registration->GetMetric());
        registration->SetOptimizer(rgsd);
        this->m_Framework->SetOptimizer(registration->GetOptimizer());
        this->m_Framework->ReAddOptimizationObservers();
  
        if (this->m_Verbose)
          std::cout << " --> SUCCESS" << std::endl;
      } // else: expected to be processed by super-class
    }
    // this->m_OptimizerConfigs: processed by superclass ...

    Superclass::ReConfigureRegistrationComponents();

    WriteFixedImageMask();
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommandWithSRC<TFramework>
::WriteFinalPDFImage(RegistrationPointer registration,
  typename BaseOptimizerType::ParametersType finalPars)
{
  typedef StochasticRankCorrelationImageToImageMetric<
    typename FrameworkType::InternalImageType,
    typename FrameworkType::InternalImageType,
    typename FrameworkType::RankImageType> SRCMetricType;
  typedef typename SRCMetricType::Pointer SRCMetricPointer;

  SRCMetricPointer src = SRCMetricType::New();
  if (registration->GetMetric()->GetNameOfClass() == src->GetNameOfClass())
  {
    if (registration && this->m_ImageOutput &&
        (this->m_ImageFileNames.size() >= 5) &&
        (this->m_ImageFileNames[4].length() > 0))
    {
      if (this->m_Verbose)
        std::cout << std::endl << "Writing final sample distribution" <<
          " to disc (level=" << registration->GetCurrentLevel() << ")." <<
          std::endl << std::endl;

      FrameworkPointer fw = this->m_Framework;

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
      sprintf(fn, this->m_ImageFileNames[4].c_str(),
        registration->GetCurrentLevel());
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
    Superclass::WriteFinalPDFImage(registration, finalPars); // delegate
  }
}


}


#endif /* ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_TXX_ */
