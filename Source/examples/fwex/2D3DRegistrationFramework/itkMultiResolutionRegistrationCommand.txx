

#ifndef ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_TXX_
#define ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_TXX_


#include <itkMultiResolutionRegistrationCommand.h>
#include <itkImageFileWriter.h>
#include <itkExtractImageFilter.h>

// support for concrete casting:
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkExhaustiveOptimizerComplete.hxx>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>


namespace itk
{


template <typename TFramework >
MultiResolutionRegistrationCommand<TFramework >
::MultiResolutionRegistrationCommand()
{
  m_Verbose = false;
  m_Logging = false;
  m_LogFileName = "";
  m_ImageOutput = false;
  m_ImageFileNames.clear();
  m_MetricConfigs.clear();
  m_MetricSels.clear();
  m_OptimizerConfigs.clear();
  m_OptimizerSels.clear();
  m_InterpolatorConfigs.clear();
}

template <typename TFramework>
MultiResolutionRegistrationCommand<TFramework>
::~MultiResolutionRegistrationCommand()
{
  ;
}

template <typename TFramework>
void MultiResolutionRegistrationCommand<TFramework>
::Execute(const Object *, const EventObject &)
{
  return;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::Execute(Object *object, const EventObject &event)
{
  if (!m_Framework || !object)
    return;

  // accept level change events only!
  if(typeid(event) != typeid(IterationEvent))
  {
    return;
  }

  // cast the basic components:
  RegistrationPointer registration = SmartPointer<RegistrationType>(
    dynamic_cast<RegistrationType *>(object));


  if (registration->GetCurrentLevel() == 0) // init stamp: time measurement
  {    
    m_Framework->SetLastTimeStamp(m_Framework->GetClock()->GetTimeStamp());
  }
  
  // current image pyramid:
  VerboseSpacingProperties(registration);

  // adjust projector's properties:
  m_Framework->AdjustProjectorOutputImageProperties();

  // image-output:
  WriteFixedImage(registration);
  WriteInputVolume(registration);
  WriteMovingImage(registration);  


  // set the final transformation parameters of previous level as initial
  // transform parameters for next level (levels > 0):
  if (registration->GetCurrentLevel() > 0)
  {
    typename BaseOptimizerType::ParametersType finalPars;
    double bestValue = 0.;
    unsigned int numberOfIterations = -1;

    // for e.g. exhaustive optimizer the init pos of next level must be
    // set manually:
    TakeOverInitialPositionForNextLevel(registration);

    // get the results (optimizer/metric-dependent):
    if (GetCurrentResults(registration, finalPars, bestValue,
        numberOfIterations))
    {
      // verbose and log summarized registration results:
      VerboseAndLogRegistrationResultAtLevel(registration,
        numberOfIterations, finalPars, bestValue);


      // write the final images for current level:
      WriteFinalMovingImage(registration, finalPars);
      WriteFinalPDFImage(registration, finalPars);


      // take over the final parameters:
      registration->GetTransform()->SetParameters(finalPars);
    }
    
    // allow re-configuration of registration components at level-changes:
    ReConfigureRegistrationComponents();
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::WriteFixedImage(RegistrationPointer registration)
{
  if (registration && m_ImageOutput)
  {
    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    char fn[2048];

    int filelevel = registration->GetCurrentLevel() + 1; // +1

    if (m_Verbose)
      std::cout << std::endl << "Writing fixed image(s) " <<
        "for level " << filelevel << std::endl << std::endl;

    // - fixed image:
    if ((m_ImageFileNames.size() >= 1) &&
        (m_ImageFileNames[0].length() > 0))
    {
      sprintf(fn, m_ImageFileNames[0].c_str(), filelevel);
      writer->SetFileName(fn);

      // NOTE: the image pyramid filter must be used because registration's
      // Initialize()-method is called after this IterationEvent (command)!
      writer->SetInput(registration->GetFixedImagePyramid()->
        GetOutput(registration->GetCurrentLevel()));
      try
      {
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing fixed image (" <<
          fn << ")." << std::endl;
      }
    }
    // - cropped fixed image:
    if ((m_ImageFileNames.size() >= 6) &&
        (m_ImageFileNames[5].length() > 0))
    {
      sprintf(fn, m_ImageFileNames[5].c_str(), filelevel);
      writer->SetFileName(fn);

      InternalImagePointer img = registration->GetFixedImagePyramid()->
        GetOutput(registration->GetCurrentLevel());

      typedef ExtractImageFilter<InternalImageType, InternalImageType> CropType;
      typedef typename CropType::Pointer CropPointer;
      typedef typename InternalImageType::RegionType CropRegionType;

      CropPointer cropper = CropType::New();

      cropper->SetInput(img);

      InternalSpacingType fixedSpacing;
      InternalRegionType fixedRegion;
      InternalSpacingType movingSpacing;
      InternalRegionType imgRegion;
      
      imgRegion = img->GetLargestPossibleRegion();
      
      m_Framework->GetCurrentImageProperties(fixedSpacing, fixedRegion,
        movingSpacing);      
            
      cropper->SetExtractionRegion(fixedRegion); // in fixedSpacing-size!
      
      try
      {        
        cropper->Update();        
        writer->SetInput(cropper->GetOutput());
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing cropped fixed image (" <<
          fn << ";" << e << ")." << std::endl;
      }
    }
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::WriteInputVolume(RegistrationPointer registration)
{
  if (registration && m_ImageOutput)
  {
    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    char fn[2048];

    int filelevel = registration->GetCurrentLevel() + 1; // +1

    if (m_Verbose)
      std::cout << std::endl << "Writing input volume " <<
        "for level " << filelevel << std::endl << std::endl;

    // - input volume image:
    if ((m_ImageFileNames.size() >= 2) &&
        (m_ImageFileNames[1].length() > 0))
    {
      sprintf(fn, m_ImageFileNames[1].c_str(), filelevel);
      writer->SetFileName(fn);

      // NOTE: the image pyramid filter must be used because registration's
      // Initialize()-method is called after this IterationEvent!
      writer->SetInput(registration->GetMovingImagePyramid()->
        GetOutput(registration->GetCurrentLevel()));
      try
      {
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing input volume image (" <<
          fn << ")." << std::endl;
      }
    }
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::WriteMovingImage(RegistrationPointer registration)
{
  if (registration && m_ImageOutput)
  {
    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    char fn[2048];

    int filelevel = registration->GetCurrentLevel() + 1; // +1

    if (m_Verbose)
      std::cout << std::endl << "Writing moving image " <<
        "for level " << filelevel << std::endl << std::endl;

    // - initial moving image:
    if ((m_ImageFileNames.size() >= 3) &&
        (m_ImageFileNames[2].length() > 0))
    {
      sprintf(fn, m_ImageFileNames[2].c_str(), filelevel);
      writer->SetFileName(fn);

      typedef typename FrameworkType::BaseTransformType TransformType;
      typedef typename TransformType::ParametersType ParametersType;

      // to get the initial moving image for current level, the transformation
      // must be set to the initial setting for current level (this is because
      // Initialize() has not been called at this point):
      ParametersType parsStore = m_Framework->GetTransformation()->
        GetParameters();

      if (registration->GetCurrentLevel() > 0)
        // not first iteration: prepared transform parameters for next level
        m_Framework->GetTransformation()->SetParameters(
          registration->GetInitialTransformParametersOfNextLevel());
      else
        // first iteration: initial transform parameters
        m_Framework->GetTransformation()->SetParameters(
          registration->GetInitialTransformParameters());      

      writer->SetInput(m_Framework->GetCurrentProjectionImage());

      try
      {
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing initial moving image (" <<
          fn << ")." << std::endl;
      }

      // NOTE: it is not easy to implement a 'initial PDF image' because
      // the multi-resolution registration internally calls ::Initialize()
      // after IterationEvent and therefore the metric has set the wrong (or
      // no) components.

      // restore parameters:
      m_Framework->GetTransformation()->SetParameters(parsStore);
    }
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::WriteFinalMovingImage(RegistrationPointer registration,
  typename BaseOptimizerType::ParametersType finalPars)
{
  if (registration && m_ImageOutput && (m_ImageFileNames.size() >= 4) &&
          (m_ImageFileNames[3].length() > 0))
  {
    if (m_Verbose)
      std::cout << std::endl << "Writing final moving image to disc " <<
        "(level=" << registration->GetCurrentLevel() << ")." <<
        std::endl << std::endl;

    // to get the final moving/PDF image for current level, the transformation
    // must be set to the result transformation parameters:
    typedef typename FrameworkType::BaseTransformType::ParametersType
      ParametersType;

    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    ParametersType parsStore = m_Framework->GetTransformation()->
      GetParameters();

    m_Framework->GetTransformation()->SetParameters(finalPars);

    // - final moving image:
    char fn[2048];
    sprintf(fn, m_ImageFileNames[3].c_str(),
      registration->GetCurrentLevel());
    writer->SetFileName(fn);

    typename ImageType::ConstPointer mov = m_Framework->
      GetCurrentProjectionImage();

    if (mov)
    {
      writer->SetInput(mov);
      try
      {
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing final moving image (" <<
          fn << "): " << e << std::endl;
      }
    }

    // restore:
    m_Framework->GetTransformation()->SetParameters(parsStore);
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::WriteFinalPDFImage(RegistrationPointer registration,
  typename BaseOptimizerType::ParametersType finalPars)
{
  if (registration && m_ImageOutput && (m_ImageFileNames.size() >= 5) &&
      (m_ImageFileNames[4].length() > 0))
  {
    std::cout << std::endl << "Writing final PDF image to disc " <<
      "(level=" << registration->GetCurrentLevel() << ")." <<
      std::endl << std::endl;

    // to get the final moving/PDF image for current level, the transformation
    // must be set to the result transformation parameters:
    typedef typename FrameworkType::BaseTransformType::ParametersType
      ParametersType;

    typedef typename FrameworkType::HistImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    ParametersType parsStore = m_Framework->GetTransformation()->
      GetParameters();

    m_Framework->GetTransformation()->SetParameters(finalPars);

    // - final PDF image:
    char fn[2048];
    sprintf(fn, m_ImageFileNames[4].c_str(),
      registration->GetCurrentLevel());
    writer->SetFileName(fn);

    typename ImageType::Pointer pdf = m_Framework->
      GetCurrentMetricHistogram(true);

    if (pdf)
    {
      writer->SetInput(pdf);
      try
      {
        writer->Update();
      }
      catch (ExceptionObject &e)
      {
        std::cerr << "Error during writing final PDF image (" <<
          fn << "): " << e << std::endl;
      }
    }

    // restore:
    m_Framework->GetTransformation()->SetParameters(parsStore);
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::VerboseSpacingProperties(RegistrationPointer registration)
{
  if (m_Verbose && registration)
  {
    std::cout << "Current multi-resolution registration level = " <<
      registration->GetCurrentLevel() << std::endl;

    typename FrameworkType::InternalSpacingType fspac;
    typename FrameworkType::InternalRegionType freg;
    typename FrameworkType::InternalSpacingType mspac;

    m_Framework->GetCurrentImageProperties(fspac, freg, mspac);

    std::cout << " (*) fixed image spacing:  " << fspac << std::endl;
    std::cout << " (*) fixed image region:";
    std::cout << "  (i) index: " << "  " << freg.GetIndex();
    std::cout << "  (ii) size: " << freg.GetSize() << std::endl;
    std::cout << " (*) moving image spacing:  " << mspac << std::endl;
    std::cout << std::endl;
  }
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::VerboseAndLogRegistrationResultAtLevel(RegistrationPointer registration,
  unsigned int numberOfIterations,
  typename BaseOptimizerType::ParametersType finalPars,
  double bestValue)
{
  std::stringstream os;
  os.str("");
  os << std::endl << std::endl;
  os << "Final Registration Parameters (level=" <<
    registration->GetCurrentLevel() << "):" << std::endl;

  os << " (*) Number of Iterations: \t\t" << numberOfIterations <<
    std::endl;
  os << " (*) Final Parameters: \t\t" << finalPars << std::endl;
  os << " (*) Best Metric: \t\t" << bestValue << std::endl;

  typename FrameworkType::RealTimeStampType levelTime =
    m_Framework->GetClock()->GetTimeStamp() - m_Framework->GetLastTimeStamp();
  os << " (*) Level Time: \t\t" << levelTime << "s" << std::endl;
  // re-new:
  m_Framework->SetLastTimeStamp(m_Framework->GetClock()->GetTimeStamp());
  os << std::endl;

  if (m_Verbose)
    std::cout << os.str() << std::endl;

  if (m_Logging && (m_LogFileName.length() > 0))
  {
    char fn[2048];

    sprintf(fn, m_LogFileName.c_str(),
      registration->GetCurrentLevel());
    m_LogFile.open(fn);

    m_LogFile << os.str();

    m_LogFile.close();
  }
}

template <typename TFramework>
bool
MultiResolutionRegistrationCommand<TFramework>
::GetCurrentResults(RegistrationPointer registration,
  typename BaseOptimizerType::ParametersType &finalPars, double &bestValue,
  unsigned int &numberOfIterations)
{
  bool validResults = false;

  /* typedefs for concrete casting */
  typedef ExhaustiveOptimizerComplete ExhOptimizerType;
  typedef ExhOptimizerType::ConstPointer ExhOptimizerConstPointer;
  typedef OnePlusOneEvolutionaryOptimizer EvolOptimizerType;
  typedef EvolOptimizerType::ConstPointer EvolOptimizerConstPointer;
  typedef NormalizedMutualInformationHistogramImageToImageMetric<
    InternalImageType, InternalImageType> NmiMetricType;
  typedef typename NmiMetricType::ConstPointer NmiMetricConstPointer;


  // concrete casts (first dummies, then references):
  ExhOptimizerConstPointer exhOptimizer = SmartPointer<const ExhOptimizerType>(
    dynamic_cast<const ExhOptimizerType *>(
    ExhOptimizerType::New().GetPointer()));
  EvolOptimizerConstPointer evolOptimizer =
    SmartPointer<const EvolOptimizerType>(
    dynamic_cast<const EvolOptimizerType *>(
    EvolOptimizerType::New().GetPointer()));
  NmiMetricConstPointer nmiMetric = SmartPointer<const NmiMetricType>(
    dynamic_cast<const NmiMetricType *>(NmiMetricType::New().GetPointer()));

  BaseOptimizerPointer optimizer = SmartPointer<BaseOptimizerType>(
    dynamic_cast<BaseOptimizerType *>(registration->GetOptimizer()));
  BaseMetricPointer metric = SmartPointer<BaseMetricType>(
    dynamic_cast<BaseMetricType *>(registration->GetMetric()));

  // final transformation parameters and final/best metric value:
  if (m_Framework->GetOptimizer()->GetNameOfClass() ==
      exhOptimizer->GetNameOfClass())
  {
    exhOptimizer = SmartPointer<ExhOptimizerType>(
      dynamic_cast<ExhOptimizerType *>(optimizer.GetPointer()));
    numberOfIterations = exhOptimizer->GetCurrentIteration() + 1;
    if (!m_Framework->GetOptimizationMinimization())
    {
      finalPars = exhOptimizer->GetMaximumMetricValuePosition();
      bestValue = exhOptimizer->GetMaximumMetricValue();
    }
    else
    {
      finalPars = exhOptimizer->GetMinimumMetricValuePosition();
      bestValue = exhOptimizer->GetMinimumMetricValue();
    }
    validResults = true;
  }
  else if (m_Framework->GetOptimizer()->GetNameOfClass() ==
           evolOptimizer->GetNameOfClass())
  {
    evolOptimizer = SmartPointer<EvolOptimizerType>(
      dynamic_cast<EvolOptimizerType *>(optimizer.GetPointer()));
    numberOfIterations = evolOptimizer->GetCurrentIteration();
    finalPars = evolOptimizer->GetCurrentPosition();
    bestValue = evolOptimizer->GetCurrentCost();
    validResults = true;
  }

  return validResults;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::TakeOverInitialPositionForNextLevel(RegistrationPointer registration)
{
  /* typedefs for concrete casting */
  typedef ExhaustiveOptimizerComplete ExhOptimizerType;
  typedef ExhOptimizerType::ConstPointer ExhOptimizerConstPointer;

  // dummy for type verification:
  ExhOptimizerConstPointer exhOptimizer = SmartPointer<const ExhOptimizerType>(
    dynamic_cast<const ExhOptimizerType *>(
    ExhOptimizerType::New().GetPointer()));

  bool currentIsExhaustive = (m_Framework->GetOptimizer()->GetNameOfClass() ==
      exhOptimizer->GetNameOfClass()); 
  unsigned int level = registration->GetCurrentLevel();
  bool nextIsExhaustive = false;
  if ((m_OptimizerSels.size() >= level && 
       m_OptimizerSels[level - 1] == "EXH") || m_OptimizerSels.size() < level)
    nextIsExhaustive = true;  

  // NOTE: Important! For exhaustive optimization (or better: exhaustive
  // cost space investigation) on continuous level the internal set 
  // 'initial transform parameters of next level' make no sense (it is the 
  // initial 'coordinate' of first search). It is better to set the 
  // 'initial transform parameters'
  // of registration again. So the cost space around initial transform
  // parameters can be investigated at the different levels.
  if (currentIsExhaustive && nextIsExhaustive)
    registration->SetInitialTransformParametersOfNextLevel(
      registration->GetInitialTransformParameters());
}

template <typename TFramework>
void 
MultiResolutionRegistrationCommand<TFramework>
::SetMetricSels(const std::vector<std::string> &sels)
{
  m_MetricSels.clear();
  for (unsigned int i = 0; i < sels.size(); i++)
    m_MetricSels.push_back(sels[i]);
}

template <typename TFramework>
const std::vector<std::string> &
MultiResolutionRegistrationCommand<TFramework>
::GetMetricSels() const
{
  return m_MetricSels;
}

template <typename TFramework>
void 
MultiResolutionRegistrationCommand<TFramework>
::SetMetricConfigs(const std::vector<std::string> &configs)
{
  m_MetricConfigs.clear();
  for (unsigned int i = 0; i < configs.size(); i++)
    m_MetricConfigs.push_back(configs[i]);
}

template <typename TFramework>
const std::vector<std::string> &
MultiResolutionRegistrationCommand<TFramework>
::GetMetricConfigs() const
{
  return m_MetricConfigs;
}

template <typename TFramework>
void 
MultiResolutionRegistrationCommand<TFramework>
::SetOptimizerSels(const std::vector<std::string> &sels)
{
  m_OptimizerSels.clear();
  for (unsigned int i = 0; i < sels.size(); i++)
    m_OptimizerSels.push_back(sels[i]);
}

template <typename TFramework>
const std::vector<std::string> &
MultiResolutionRegistrationCommand<TFramework>
::GetOptimizerSels() const
{
  return m_OptimizerSels;
}

template <typename TFramework>
void 
MultiResolutionRegistrationCommand<TFramework>
::SetOptimizerConfigs(const std::vector<std::string> &configs)
{
  m_OptimizerConfigs.clear();
  for (unsigned int i = 0; i < configs.size(); i++)
    m_OptimizerConfigs.push_back(configs[i]);
}

template <typename TFramework>
const std::vector<std::string> &
MultiResolutionRegistrationCommand<TFramework>
::GetOptimizerConfigs() const
{
  return m_OptimizerConfigs;
}

template <typename TFramework>
void 
MultiResolutionRegistrationCommand<TFramework>
::SetInterpolatorConfigs(const std::vector<std::string> &configs)
{
  m_InterpolatorConfigs.clear();
  for (unsigned int i = 0; i < configs.size(); i++)
    m_InterpolatorConfigs.push_back(configs[i]);
}

template <typename TFramework>
const std::vector<std::string> &
MultiResolutionRegistrationCommand<TFramework>
::GetInterpolatorConfigs() const
{
  return m_InterpolatorConfigs;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommand<TFramework>
::ReConfigureRegistrationComponents()
{  
  RegistrationPointer registration = m_Framework->GetRegistration();
  unsigned int level = registration->GetCurrentLevel(); 
  if (level < 1)
    return;
  
  bool succ;
  
  // metric
  if (m_MetricSels.size() >= level && 
      m_MetricConfigs.size() >= m_MetricSels.size()) // need config as well!
  {        
    if (m_MetricSels[level - 1] == "NMI")
    {
      if (m_Verbose)
        std::cout << std::endl << "Resetting metric component for level " <<
          (level + 1) << ": " << m_MetricSels[level - 1];
          
      typedef NormalizedMutualInformationHistogramImageToImageMetric<
        InternalImageType, InternalImageType> NmiMetricType;
      typename NmiMetricType::Pointer nmi = NmiMetricType::New();
      nmi->SetFixedImage(registration->GetFixedImage());
      nmi->SetMovingImage(registration->GetMovingImage());
      nmi->SetTransform(registration->GetTransform());
      nmi->SetInterpolator(registration->GetInterpolator());
      typename InternalImageType::RegionType freg = registration->GetMetric()->
        GetFixedImageRegion();
      nmi->SetFixedImageRegion(freg);      
      registration->SetMetric(nmi);
      m_Framework->SetMetric(registration->GetMetric());

      if (m_Verbose)
        std::cout << " --> SUCCESS" << std::endl;
    } // else: expected to be processed by sub-classes
  }
  if (m_MetricConfigs.size() >= level) // sufficient info available
  {    
    if (m_Verbose)
      std::cout << std::endl << "Reconfiguring metric component for level " <<
        (level + 1) << ": " << m_MetricConfigs[level - 1];
    succ = m_Framework->ConfigureMetric(m_MetricConfigs[level - 1]);
    if (m_Verbose)
    {
      if (succ)
        std::cout << " --> SUCCESS" << std::endl << std::endl;
      else
        std::cout << " --> FAILURE" << std::endl << std::endl;
    }
  }
  
  // optimizer
  if (m_OptimizerSels.size() >= level && 
      m_OptimizerConfigs.size() >= m_OptimizerSels.size()) // need config, too!
  {
    if (m_Verbose)
      std::cout << std::endl << "Resetting optimizer component for level " <<
        (level + 1) << ": " << m_OptimizerSels[level - 1];    
    if (m_OptimizerSels[level - 1] == "EXH")
    {
      typedef ExhaustiveOptimizerComplete ExhOptimizerType;
      ExhOptimizerType::Pointer exh = ExhOptimizerType::New();
      exh->SetCostFunction(registration->GetMetric());      
      registration->SetOptimizer(exh);
      m_Framework->SetOptimizer(registration->GetOptimizer());
      m_Framework->ReAddOptimizationObservers();
      if (m_Verbose)
        std::cout << " --> SUCCESS" << std::endl;
    }
    else if (m_OptimizerSels[level - 1] == "EVOL")
    {
      typedef OnePlusOneEvolutionaryOptimizer EvolOptimizerType;
      EvolOptimizerType::Pointer evol = EvolOptimizerType::New();
      evol->SetCostFunction(registration->GetMetric());
      registration->SetOptimizer(evol);
      m_Framework->SetOptimizer(registration->GetOptimizer());
      m_Framework->ReAddOptimizationObservers();
      if (m_Verbose)
        std::cout << " --> SUCCESS" << std::endl;
    }
  }
  if (m_OptimizerConfigs.size() >= level) // sufficient info available
  {    
    if (m_Verbose)
      std::cout << std::endl << "Reconfiguring optimizer component for level "<<
        (level + 1) << ": " << m_OptimizerConfigs[level - 1];
    succ = m_Framework->ConfigureOptimizer(m_OptimizerConfigs[level - 1]);
    if (m_Verbose)
    {
      if (succ)
        std::cout << " --> SUCCESS" << std::endl << std::endl;
      else
        std::cout << " --> FAILURE" << std::endl << std::endl;
    }
  }
  
  // interpolator
  if (m_InterpolatorConfigs.size() >= level) // sufficient info available
  {    
    if (m_Verbose)
      std::cout << std::endl << "Reconfiguring interpolator component for " <<
        "level " << (level + 1) << ": " << m_InterpolatorConfigs[level - 1];
    succ = m_Framework->ConfigureInterpolator(m_InterpolatorConfigs[level - 1]);
    if (m_Verbose)
    {
      if (succ)
        std::cout << " --> SUCCESS" << std::endl << std::endl;
      else
        std::cout << " --> FAILURE" << std::endl << std::endl;
    }
  }
}


}


#endif /* ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_TXX_ */
