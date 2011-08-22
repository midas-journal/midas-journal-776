

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_TXX_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_TXX_


#include <itkMultiResolutionImage2D3DRegistrationMethod.h>

#include <vector>

#include <itkBSplineInterpolateImageFunction.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkHistogramToEntropyImageFilter.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>

// to support the configure-methods:
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkRayCastPerspectiveProjectionImageFilter.h>
#include <itkSplatPerspectiveProjectionImageFilter.h>
#include <itkEuler3DTransform.h>
#include <itkAffineTransform.h>
#include <itkNormalVariateGenerator.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkExhaustiveOptimizerComplete.hxx>


namespace itk
{


template <typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::MultiResolutionImage2D3DRegistrationMethod()
{

  m_FixedImage = NULL;
  m_RescaledFixedImage = NULL;
  m_FixedImageRegion.SetIndex(0, 0);
  m_FixedImageRegion.SetIndex(1, 0);
  m_FixedImageRegion.SetIndex(2, 0);
  m_FixedImageRegion.SetSize(0, 0);
  m_FixedImageRegion.SetSize(1, 0);
  m_FixedImageRegion.SetSize(2, 1);
  m_MovingImage = NULL;
  m_XrayImageRegion.SetIndex(0, 0);
  m_XrayImageRegion.SetIndex(1, 0);
  m_XrayImageRegion.SetIndex(2, 0);
  m_XrayImageRegion.SetSize(0, 0);
  m_XrayImageRegion.SetSize(1, 0);
  m_XrayImageRegion.SetSize(2, 1);
  
  m_FixedRescalingActivated = false;
  m_FixedRescalingMin = 0.;
  m_FixedRescalingMax = 255.;

  m_Registration = NULL;
  m_Transformation = NULL;
  m_Interpolator = NULL;
  m_ProjectionFilter = NULL;
  m_Metric = NULL;
  m_Optimizer = NULL;  

  this->SetNumberOfRequiredOutputs(1); // @see Initialize()
  
  m_Clock = RealTimeClock::New();
  m_LastTimeStamp = m_Clock->GetTimeStamp();
  
  m_InitialTransformationParamters.Fill(0);
  
  m_OptimizerObservers.clear();
}

template <typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::~MultiResolutionImage2D3DRegistrationMethod()
{
  ;
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  std::cout << " -- 2D/3D Registration Framework -- " << std::endl;

  std::cout << std::endl << " (*) Framework:" << std::endl;
  std::cout << "Fixed Image: " << m_FixedImage << std::endl;
  std::cout << "Moving Image: " << m_MovingImage << std::endl;
  std::cout << "X-ray Image Region: " << m_XrayImageRegion << std::endl;
  std::cout << "X-ray Image Spacing: " << m_XrayImageSpacing << std::endl;

  std::cout << std::endl << " (*) Internal Registration:" << std::endl;
  std::cout << "Internal Registration: " << m_Registration << std::endl;

  std::cout << std::endl << " (*) Transformation:" << std::endl;
  std::cout << "Transformation: " << m_Transformation << std::endl;

  std::cout << std::endl << " (*) Interpolator (Projector):" << std::endl;
  std::cout << "Interpolator: " << m_Interpolator << std::endl;
  std::cout << "Internal Projector: " << m_ProjectionFilter << std::endl;

  std::cout << std::endl << " (*) Metric:" << std::endl;
  std::cout << "Metric: " << m_Metric << std::endl;

  std::cout << std::endl << " (*) Optimizer:" << std::endl;
  std::cout << "Optimizer: " << m_Optimizer << std::endl;
  std::cout << "Optimization Minimization: " << m_OptimizationMinimization << std::endl;
}

template <typename TInternalPixelType, typename TScalarType>
const typename MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, 
  TScalarType>::BaseTransformOutputType *
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GetOutput() const
{
  return static_cast<const BaseTransformOutputType *>(
    this->ProcessObject::GetOutput(0));
}

template <typename TInternalPixelType, typename TScalarType>
DataObject::Pointer 
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::MakeOutput(unsigned int idx)
{
  switch (idx)
  {
    case 0:
      return static_cast<DataObject *>(
        BaseTransformOutputType::New().GetPointer());
      break;
      
    default:
      itkExceptionMacro(<< "MakeOutput request for an output number larger " <<
        "than the expected number of outputs");
      return 0;      
  }
}

template <typename TInternalPixelType, typename TScalarType>
unsigned long 
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GetMTime() const
{
  unsigned long mtime = Superclass::GetMTime();
  unsigned long mt;
  
  if (m_Registration) // handles most of the registration-components internally
  {
    mt = m_Registration->GetMTime();
    mtime = (mt > mtime ? mt : mtime);
  }
  
  if (m_ProjectionFilter)
  {
    mt = m_Registration->GetMTime();
    mtime = (mt > mtime ? mt : mtime);
  }
  
  return mtime;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::Initialize()
{
  if (!m_FixedImage || !m_MovingImage || !m_Metric || !m_Optimizer ||
      !m_Transformation || !m_ProjectionFilter)
    return false;

  m_Registration = RegistrationType::New();

  // tie the components together:
  m_Registration->SetFixedImage(m_FixedImage); // eventually overridden  
  m_Registration->SetFixedImageRegion(m_FixedImageRegion);
  m_Registration->SetMovingImage(m_MovingImage);

  m_Registration->SetMetric(m_Metric);

  m_Registration->SetOptimizer(m_Optimizer);

  m_Registration->SetTransform(m_Transformation);

  m_Interpolator = BaseInterpolatorType::New();
  m_ProjectionFilter->SetInput(m_MovingImage);
  m_ProjectionFilter->SetTransform(m_Transformation);
  m_Interpolator->SetPerspectiveProjectionImageFilter(m_ProjectionFilter);
  m_Registration->SetInterpolator(m_Interpolator);

  BaseTransformOutputPointer transformDecorator = 
    static_cast<BaseTransformOutputType *>(
      this->MakeOutput(0).GetPointer());
  transformDecorator->Set(m_Transformation);
  this->ProcessObject::SetNthOutput(0, transformDecorator.GetPointer());

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::Deinitialize()
{
  if (!m_Registration)
    return;
  
  if (m_Registration->GetTransform())
    m_Registration->GetTransform()->RemoveAllObservers();
  m_Registration->SetTransform(NULL);

  m_Registration->SetFixedImage(NULL);
  m_Registration->SetMovingImage(NULL); 
    
  m_Registration->SetInterpolator(NULL);
  m_Interpolator = NULL;

  m_ProjectionFilter = NULL;
  
  if (m_Registration->GetMetric())
  {
    //m_Registration->GetMetric()->SetTransform(NULL);
    m_Registration->SetMetric(NULL);
  }

  if (m_Registration->GetOptimizer())
  {
    m_Registration->GetOptimizer()->RemoveAllObservers();
    m_Registration->GetOptimizer()->SetCostFunction(NULL);
    m_Registration->SetOptimizer(NULL);
  }

  m_Registration->RemoveAllObservers();

  m_Registration = NULL;

  this->RemoveAllObservers();
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::CallArgsInitialize(std::string s, itksys::CommandLineArguments *cla)
{
  // find parts and store in string-vector:
  std::vector<std::string> parts;
  std::string::size_type off = static_cast<std::string::size_type>(0);
  std::string::size_type pos = std::string::npos;
  std::string part;
  do
  {
    pos = s.find(" ", off);
    if (pos != std::string::npos)
    {
      part = s.substr(off, pos - off);
      off = pos + static_cast<std::string::size_type>(1);
    }
    else
    {
      part = s.substr(off);
    }
    if (part.length() > 0)
      parts.push_back(part);
  } while (pos != std::string::npos);

  unsigned int argc = parts.size();
  std::vector<char *> argv;

  argv.push_back(const_cast<char *>("dummy")); // 'binary name'
  for (unsigned int i = 0; i < argc; ++i)
    argv.push_back(const_cast<char *>(parts[i].c_str()));
  argv.push_back(0);

  cla->StoreUnusedArguments(true);
  cla->Initialize(++argc, &argv[0]);
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureImages(std::string config)
{

  if (!m_Registration || !m_MovingImage || !m_FixedImage)
    return false;

  itksys::CommandLineArguments cla;
  itksys_stl::vector<double> argVolumeOrigin;
  InternalPointType volumeorigin;
  itksys_stl::vector<double> argFixedRescaling;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("volumeorigin", A_M_SPAC, &argVolumeOrigin, "");
  cla.AddArgument("fixedrescaling", A_M_SPAC, &argFixedRescaling, ""); 
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ... 

  // process found arguments:
  if (argVolumeOrigin.size() == 3)
  {
    for (unsigned int i = 0; i < 3; ++i)
      volumeorigin[i] = argVolumeOrigin[i];
  }

  if (argFixedRescaling.size() >= 2)
  {
    m_FixedRescalingMin = argFixedRescaling[0];
    m_FixedRescalingMax = argFixedRescaling[1];
    m_FixedRescalingActivated = true; 
  }
  else
    m_FixedRescalingActivated = false;

  // set up images:
  m_MovingImage->SetOrigin(volumeorigin);

  // rescale on demand and set it in framework:
  if (m_FixedRescalingActivated)
  {
    RescaleFixedImageInRegion(m_FixedRescalingMin, m_FixedRescalingMax);
    m_Registration->SetFixedImage(m_RescaledFixedImage);
  }

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureImagePyramids(const std::string config)
{
  if (!m_Registration)
    return false;

  itksys::CommandLineArguments cla;
  int pyramidLevels = -1;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("levels", A_SPAC, &pyramidLevels, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // process found arguments:
  if (pyramidLevels < 1)
    return false; // not specified

  typename RegistrationType::ScheduleType fschedule(pyramidLevels, 3);
  fschedule.Fill(1);
  typename RegistrationType::ScheduleType mschedule(pyramidLevels, 3);
  mschedule.Fill(1);
  for (int l = 0; l < pyramidLevels; l++)
  {
    unsigned int fact = static_cast<unsigned int>(pow((double)2, (int)(pyramidLevels - l - 1)));
    for (unsigned int d = 0; d < 3; d++)
    {
      mschedule[l][d] = fact;
      if (d < 2) // do not scale fixed image's virtual 3rd dimension!
        fschedule[l][d] = fact;
    }
  }
  m_Registration->SetSchedules(fschedule, mschedule);  

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureTransformation(const std::string config)
{
  if (!m_Transformation)
    return false;

  /* typedefs for concrete casting */
  typedef Euler3DTransform<ScalarType> EulerTransformType;
  typedef typename EulerTransformType::Pointer EulerTransformPointer;
  typedef AffineTransform<ScalarType> AffineTransformType;
  typedef typename AffineTransformType::Pointer AffineTransformPointer;

  itksys::CommandLineArguments cla;
  typename BaseTransformType::InputPointType cor;
  itksys_stl::vector<double> argCor;
  itksys_stl::vector<double> rot;
  typename BaseTransformType::OutputVectorType trans;
  itksys_stl::vector<double> argTrans;
  itksys_stl::vector<double> argScale;
  typename BaseTransformType::OutputVectorType scale;

  const double degToRad = acos((double)-1) / 180.;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("cor", A_M_SPAC, &argCor, "");
  cla.AddArgument("rot", A_M_SPAC, &rot, "");
  cla.AddArgument("transl", A_M_SPAC, &argTrans, "");
  cla.AddArgument("scale", A_M_SPAC, &argScale, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // dummies for type verification:
  EulerTransformPointer eulerDummy = EulerTransformType::New();
  AffineTransformPointer affineDummy = AffineTransformType::New();

  // process found arguments:
  if (m_Transformation->GetNameOfClass() == eulerDummy->GetNameOfClass() ||
      m_Transformation->GetNameOfClass() == affineDummy->GetNameOfClass())
  {
    if (argCor.size() == 3)
    {
      for (unsigned int i = 0; i < 3; ++i)
        cor[i] = argCor[i];
    }
    else
      return false; // wrong number of arguments
    if (rot.size() == 3)
    {
      for (unsigned int i = 0; i < 3; ++i)
        rot[i] = rot[i] * degToRad;
    }
    else
      return false; // wrong number of arguments
    if (argTrans.size() == 3)
    {
      for (unsigned int i = 0; i < 3; ++i)
        trans[i] = argTrans[i];
    }
    else
      return false; // wrong number of arguments

    if (m_Transformation->GetNameOfClass() == affineDummy->GetNameOfClass())
    {
      if (argScale.size() == 3)
      {
        for (unsigned int i = 0; i < 3; ++i)
          scale[i] = argScale[i];
      }
      else if (argScale.size() != 0)
        return false; // wrong number of arguments
    }
  }

  // set-up transformation:
  if (m_Transformation->GetNameOfClass() == eulerDummy->GetNameOfClass())
  {
    EulerTransformPointer transform = SmartPointer<EulerTransformType>(
      static_cast<EulerTransformType *>(m_Transformation.GetPointer()));

    transform->SetComputeZYX(true);

    if (m_MovingImage) // add origin to center of rotation
    { InternalPointType origin = m_MovingImage->GetOrigin();

      for (int z = 0; z < origin.PointDimension; z++)
        cor[z] = cor[z] + origin[z];
    }
    transform->SetCenter(cor);

    transform->SetRotation(rot[0], rot[1], rot[2]);
    transform->SetTranslation(trans);
  }
  else if (m_Transformation->GetNameOfClass() == affineDummy->GetNameOfClass())
  {
    AffineTransformPointer transform = SmartPointer<AffineTransformType>(
      static_cast<AffineTransformType *>(m_Transformation.GetPointer()));
    typename AffineTransformType::OutputVectorType rotAxis;

    if (m_MovingImage) // add origin to center of rotation
    { InternalPointType origin = m_MovingImage->GetOrigin();

      for (int z = 0; z < origin.PointDimension; z++)
        cor[z] = cor[z] + origin[z];
    }
    transform->SetCenter(cor);

    if (argScale.size() == 3)
      transform->Scale(scale, false);

    // z-axis:
    rotAxis[0] = 0;
    rotAxis[1] = 0;
    rotAxis[2] = 1;
    transform->Rotate3D(rotAxis, rot[2], true);
    // y-axis:
    rotAxis[0] = 0;
    rotAxis[1] = 1;
    rotAxis[2] = 0;
    transform->Rotate3D(rotAxis, rot[1], true);
    // x-axis:
    rotAxis[0] = 1;
    rotAxis[1] = 0;
    rotAxis[2] = 0;
    transform->Rotate3D(rotAxis, rot[0], true);

    transform->SetTranslation(trans);
  }

  m_Registration->SetInitialTransformParameters(m_Transformation->
    GetParameters());

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureOptimizer(const std::string config)
{
  if (!m_Optimizer || !m_Transformation)
    return false;

  /* typedefs for concrete casting */
  typedef ExhaustiveOptimizerComplete ExhOptimizerType;
  typedef ExhOptimizerType::Pointer ExhOptimizerPointer;
  typedef OnePlusOneEvolutionaryOptimizer EvolOptimizerType;
  typedef EvolOptimizerType::Pointer EvolOptimizerPointer;
  typedef Statistics::NormalVariateGenerator NVGeneratorType;
  typedef NVGeneratorType::Pointer NVGeneratorPointer;

  itksys::CommandLineArguments cla;
  int min = -1;
  Optimizer::ScalesType weights;
  itksys_stl::vector<double> argScales;
  double steplen = -1;
  ExhOptimizerType::StepsType stepcount;
  itksys_stl::vector<int> argStepCount;
  int iterations = 0;
  int nvgseed = 12345;
  double radius = 1.;
  double gfact = -1;
  double sfact = -1;
  double epsilon = 1e-4;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("min", A_SPAC, &min, "");
  cla.AddArgument("scales", A_M_SPAC, &argScales, "");
  cla.AddArgument("steplen", A_SPAC, &steplen, "");
  cla.AddArgument("stepcount", A_M_SPAC, &argStepCount, "");
  cla.AddArgument("iterations", A_SPAC, &iterations, "");
  cla.AddArgument("nvgseed", A_SPAC, &nvgseed, "");
  cla.AddArgument("radius", A_SPAC, &radius, "");
  cla.AddArgument("gfact", A_SPAC, &gfact, "");
  cla.AddArgument("sfact", A_SPAC, &sfact, "");
  cla.AddArgument("epsilon", A_SPAC, &epsilon, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // dummies for type verification:
  ExhOptimizerPointer exhDummy = ExhOptimizerType::New();
  EvolOptimizerPointer evolDummy = EvolOptimizerType::New();

  // process found arguments:
  if ((min != 0) && (min != 1))
    return false; // wrong value
  if (argScales.size() == m_Transformation->GetNumberOfParameters())
  {
    weights.SetSize(argScales.size());
    for (unsigned int i = 0; i < argScales.size(); ++i)
      weights[i] = argScales[i];
  }
  else
    return false; // wrong number of arguments
  if (m_Optimizer->GetNameOfClass() == exhDummy->GetNameOfClass())
  {
    if (steplen < 0)
      return false; // wrong value
    if (argStepCount.size() == argScales.size())
    {
      stepcount.SetSize(argStepCount.size());
      for (unsigned int i = 0; i < argStepCount.size(); ++i)
        stepcount[i] = argStepCount[i];
    }
    else
      return false; // wrong number of arguments
  }
  else if (m_Optimizer->GetNameOfClass() == evolDummy->GetNameOfClass())
  {
    if (iterations <= 0)
      return false; // not allowed
    if (radius <= 0.0)
      return false; // not allowed
  }


  // set-up optimizer:
  // - common for all optimizers:
  m_OptimizationMinimization = min;
  m_Optimizer->SetScales(weights);

  if (m_Optimizer->GetNameOfClass() == exhDummy->GetNameOfClass())
  {
    ExhOptimizerPointer optimizer;

    optimizer = SmartPointer<ExhOptimizerType>(static_cast<ExhOptimizerType *>(
      m_Optimizer.GetPointer()));

    // min: not used for optimizer, but internally stored in framework;
    optimizer->SetStepLength(steplen);
    optimizer->SetNumberOfSteps(stepcount);
  }
  else if (m_Optimizer->GetNameOfClass() == evolDummy->GetNameOfClass())
  {
    EvolOptimizerPointer optimizer;

    optimizer = SmartPointer<EvolOptimizerType>(
      static_cast<EvolOptimizerType *>(m_Optimizer.GetPointer()));

    NVGeneratorPointer gen = NVGeneratorType::New();

    gen->Initialize(nvgseed);
    optimizer->SetNormalVariateGenerator(gen);

    optimizer->SetMaximize(!m_OptimizationMinimization);

    optimizer->SetMaximumIteration(iterations);

    // the real arguments:
    optimizer->Initialize(radius, gfact, sfact);
    optimizer->SetEpsilon(epsilon);
  }

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureMetric(const std::string config)
{
  if (!m_Metric)
    return false;

  /* typedefs for concrete casting */
  typedef NormalizedMutualInformationHistogramImageToImageMetric<
    InternalImageType, InternalImageType> NmiMetricType;
  typedef typename NmiMetricType::Pointer NmiMetricPointer;

  itksys::CommandLineArguments cla;
  typename NmiMetricType::HistogramType::SizeType nmibins;
  typename NmiMetricType::MeasurementVectorType nmilb;
  typename NmiMetricType::MeasurementVectorType nmiub;
  itksys_stl::vector<int> argNmiBins;
  itksys_stl::vector<double> argNmiLb;
  itksys_stl::vector<double> argNmiUb;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("histbins", A_M_SPAC, &argNmiBins, "");
  cla.AddArgument("lowerbounds", A_M_SPAC, &argNmiLb, "");
  cla.AddArgument("upperbounds", A_M_SPAC, &argNmiUb, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // dummies for type verification:
  NmiMetricPointer nmiDummy = NmiMetricType::New();

  // process found arguments:
  if (m_Metric->GetNameOfClass() == nmiDummy->GetNameOfClass())
  {
    if (argNmiBins.size() == 2)
    {
      for (unsigned int i = 0; i < argNmiBins.size(); ++i)
        nmibins[i] = argNmiBins[i];
    }
    else
      return false; // wrong number of arguments
    if (argNmiLb.size() == 2)
    {
      for (unsigned int i = 0; i < argNmiLb.size(); ++i)
        nmilb[i] = argNmiLb[i];
    }
    else
      return false; // wrong number of arguments
    if (argNmiUb.size() == 2)
    {
      for (unsigned int i = 0; i < argNmiUb.size(); ++i)
        nmiub[i] = argNmiUb[i];
    }
    else
      return false; // wrong number of arguments
  }


  // set-up metric:
  if (m_Metric->GetNameOfClass() == nmiDummy->GetNameOfClass())
  {
    NmiMetricPointer metric;

    metric = SmartPointer<NmiMetricType>(static_cast<NmiMetricType *>(
      m_Metric.GetPointer()));

    metric->SetHistogramSize(nmibins);
    metric->SetLowerBound(nmilb);
    metric->SetUpperBound(nmiub);
  }

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ConfigureInterpolator(const std::string config)
{
  if (!m_Interpolator || !m_ProjectionFilter)
    return false;

  /* typedefs for concrete casting */
  typedef double CoordinateRepresentationType;
  typedef Collector::SumCollector<InternalPixelType> CollectorType;
  typedef RayCastPerspectiveProjectionImageFilter<InternalImageType,
    InternalImageType, CoordinateRepresentationType, CollectorType>
    RayCastProjectionFilterType;
  typedef typename RayCastProjectionFilterType::Pointer
    RayCastProjectionFilterPointer;
  typedef SplatPerspectiveProjectionImageFilter<InternalImageType,
    InternalImageType, CoordinateRepresentationType, CollectorType>
    SplatProjectionFilterType;
  typedef typename SplatProjectionFilterType::Pointer
    SplatProjectionFilterPointer;

  itksys::CommandLineArguments cla;
  int jobs = 1;
  double thresh = 0;
  typename BaseProjectionFilterType::IntensityTransferFunctionNodesType transfer;
  itksys_stl::vector<double> argTransfer;
  typename BaseProjectionFilterType::IntensityInterpolationEnumType transferint;
  int argTransferInt = 0; // BaseProjectionFilterType::II_NONE
  typename BaseProjectionFilterType::InputPointType focalpoint;
  itksys_stl::vector<double> argFocalPoint;
  itksys_stl::vector<double> argRescale;
  typename SplatProjectionFilterType::FocalOffsetComputationEnumType
    splatfoffcomp = SplatProjectionFilterType::FOC_NONE;
  int argSplatOffComp = 0; // SplatProjectionFilterType::FOC_NONE
  double splatfoffparsigma = 1.;
  itksys_stl::vector<double> argFoffPar;
  int splatprecompvals = 50;
  typename SplatProjectionFilterType::SplatAssignmentEnumType splatpixassign =
    SplatProjectionFilterType::SA_NEAREST;
  int argPixAssign = 0; // SplatProjectionFilterType::SA_NEAREST
  bool splatfoffdeterm = false;
  double rayintstepsize = 1.;
  int rayintmethod = 0; // NearestNeighbour

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("jobs", A_SPAC, &jobs, "");
  cla.AddArgument("threshold", A_SPAC, &thresh, "");
  cla.AddArgument("transferfunc", A_M_SPAC, &argTransfer, "");
  cla.AddArgument("transferinterp", A_SPAC, &argTransferInt, "");
  cla.AddArgument("focalpoint", A_M_SPAC, &argFocalPoint, "");
  cla.AddArgument("rescale", A_M_SPAC, &argRescale, "");
  cla.AddArgument("wobbmode", A_SPAC, &argSplatOffComp, "");
  cla.AddArgument("wobbsigma", A_M_SPAC, &splatfoffparsigma, "");
  cla.AddArgument("wobbcount", A_SPAC, &splatprecompvals, "");
  cla.AddArgument("splatpixassign", A_SPAC, &argPixAssign, "");
  cla.AddArgument("wobbdeterm", A_SPAC, &splatfoffdeterm, "");
  cla.AddArgument("raystepsize", A_SPAC, &rayintstepsize, "");
  cla.AddArgument("raystepinterp", A_SPAC, &rayintmethod, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // dummies for type verification:
  RayCastProjectionFilterPointer raycDummy = RayCastProjectionFilterType::New();
  SplatProjectionFilterPointer splatDummy = SplatProjectionFilterType::New();

  // process found arguments:
  if (jobs <= 0)
    jobs = 0;
  if ((argTransfer.size() > 0) &&
      ((argTransfer.size() % 2) == 1) &&
      ((int)argTransfer[0] == (int)((argTransfer.size() - 1) / 2)))
  {
    for (unsigned int i = 1; i < argTransfer.size(); ++i)
    {
      if ((i % 2) == 0)
        transfer.insert(std::make_pair(argTransfer[i - 1],
          argTransfer[i])); // add it!
    }
  }
  else if (argTransfer.size() == 0)
    transfer.clear();
  else
    return false; // invalid transfer function
  if (argTransferInt == 0)
    transferint = BaseProjectionFilterType::II_NONE;
  else if (argTransferInt == 1)
    transferint = BaseProjectionFilterType::II_STEPWISE;
  else if (argTransferInt == 2)
    transferint = BaseProjectionFilterType::II_LINEAR;
  else
    return false; // wrong value
  if (argFocalPoint.size() == 3)
  {
    for (unsigned int i = 0; i < argFocalPoint.size(); ++i)
      focalpoint[i] = argFocalPoint[i];
  }
  else
    return false; // wrong number of arguments
  if (argRescale.size() != 0 && argRescale.size() != 2)
    return false; // wrong number of arguments (not specified allowed, optional)
  else if (argRescale.size() == 2 && argRescale[0] >= argRescale[1])
    return false; // first arg (min) must be smaller than second arg (max)
  if (m_ProjectionFilter->GetNameOfClass() == raycDummy->GetNameOfClass())
  {
    if ((rayintmethod < 0) || (rayintmethod > 3))
      return false; // wrong value
  }
  else if (m_ProjectionFilter->GetNameOfClass() == splatDummy->GetNameOfClass())
  {
    if (argSplatOffComp == 0)
      splatfoffcomp = SplatProjectionFilterType::FOC_NONE;
    else if (argSplatOffComp == 1)
      splatfoffcomp = SplatProjectionFilterType::FOC_PRE_CONTINUOUS;
    else if (argSplatOffComp == 2)
      splatfoffcomp = SplatProjectionFilterType::FOC_PRE_RANDOM;
    else if (argSplatOffComp == 3)
      splatfoffcomp = SplatProjectionFilterType::FOC_RANDOM;
    else
      return false; // wrong value
    if (splatprecompvals < 0)
      return false; // bad value
    if (argPixAssign == 0)
      splatpixassign = SplatProjectionFilterType::SA_NEAREST;
    else if (argPixAssign == 1)
      splatpixassign = SplatProjectionFilterType::SA_LINEAR;
    else
      return false; // wrong value
  }


  // set-up interpolator:
  // - common for all projection types:
  m_ProjectionFilter->SetNumberOfThreads(jobs);
  m_ProjectionFilter->SetThreshold(thresh);
  m_ProjectionFilter->SetIntensityTransferFunctionNodes(transfer);
  m_ProjectionFilter->SetIntensityInterpolation(transferint);
  m_ProjectionFilter->SetFocalPoint(focalpoint);

  typename BaseProjectionFilterType::OutputOriginType outOrigin;
  for (int z = 0; z < 3; z++)
    outOrigin[z] = m_XrayImageRegion.GetIndex()[z] * m_XrayImageSpacing[z];
  m_ProjectionFilter->SetOutputOrigin(outOrigin);
  m_ProjectionFilter->SetOutputSize(m_XrayImageRegion.GetSize());
  m_ProjectionFilter->SetOutputSpacing(m_XrayImageSpacing);

  // OPTIONAL intermediate filtering:
  if (argRescale.size() == 2)
  {
    typedef RescaleIntensityImageFilter<
      typename BaseProjectionFilterType::OutputImageType,
      typename BaseProjectionFilterType::OutputImageType> RescaleFilterType;

    typename RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();

    // pre-configuration:
    rescaleFilter->SetOutputMinimum(argRescale[0]);
    rescaleFilter->SetOutputMaximum(argRescale[1]);

    m_Interpolator->SetFindMinMaxIntensityMode(false);
    m_Interpolator->SetIntermediateFilter(rescaleFilter);
  }
  else
  {
    m_Interpolator->SetFindMinMaxIntensityMode(false);
    m_Interpolator->SetIntermediateFilter(NULL);
  }

  if (m_ProjectionFilter->GetNameOfClass() == raycDummy->GetNameOfClass())
  {
    RayCastProjectionFilterPointer raycastFilter =
      SmartPointer<RayCastProjectionFilterType>(
      static_cast<RayCastProjectionFilterType *>(
      m_ProjectionFilter.GetPointer()));

    raycastFilter->SetRayStepSize(rayintstepsize);
    if (rayintmethod == 0)
    {
      typedef itk::NearestNeighborInterpolateImageFunction<InternalImageType,
        double> InterpolatorType;

      typename InterpolatorType::Pointer raycastInterpolator =
        InterpolatorType::New();

      raycastFilter->SetInterpolator(raycastInterpolator);
    }
    if (rayintmethod == 1)
    {
      typedef itk::LinearInterpolateImageFunction<InternalImageType,
        double> InterpolatorType;

      typename InterpolatorType::Pointer raycastInterpolator =
        InterpolatorType::New();

      raycastFilter->SetInterpolator(raycastInterpolator);
    }
    if (rayintmethod == 2)
    {
      typedef itk::BSplineInterpolateImageFunction<InternalImageType,
        double> InterpolatorType;

      typename InterpolatorType::Pointer raycastInterpolator =
        InterpolatorType::New();

      raycastFilter->SetInterpolator(raycastInterpolator);
    }
    if (rayintmethod == 3)
    {
      typedef itk::Function::WelchWindowFunction<5> WindowFunctionType;
      typedef itk::ConstantBoundaryCondition<InternalImageType>
        BoundaryConditionType;
      typedef itk::WindowedSincInterpolateImageFunction<InternalImageType, 5,
        WindowFunctionType, BoundaryConditionType, double> InterpolatorType;

      typename InterpolatorType::Pointer raycastInterpolator =
        InterpolatorType::New();

      raycastFilter->SetInterpolator(raycastInterpolator);
    }
  }
  else if (m_ProjectionFilter->GetNameOfClass() == splatDummy->GetNameOfClass())
  {
    SplatProjectionFilterPointer splatFilter =
        SmartPointer<SplatProjectionFilterType>(
        static_cast<SplatProjectionFilterType *>(
        m_ProjectionFilter.GetPointer()));

    splatFilter->SetFocalOffsetComputation(splatfoffcomp);
    splatFilter->SetFocalRandomOffsetSigma(splatfoffparsigma);
    splatFilter->SetFocalRandomOffsetQuantity(splatprecompvals);
    splatFilter->SetSplatAssignment(splatpixassign);
    splatFilter->SetComputeFocalOffsetDeterministically(splatfoffdeterm);
  }

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ComputeVolumeTranslationByReferencePoint(const std::string config)
{
  if (!m_MovingImage || !m_Interpolator || !m_ProjectionFilter ||
      !m_Transformation)
    return false;

  itksys::CommandLineArguments cla;
  typename BaseProjectionFilterType::InputPointType inplanepoint;
  itksys_stl::vector<double> argInplanePoint;
  typename BaseProjectionFilterType::InputPointType referencepoint;
  itksys_stl::vector<double> argReferencePoint;
  bool applytransl = false;

  // parse configure string:
  CallArgsInitialize(config, &cla);
  cla.AddArgument("inplane", A_M_SPAC, &argInplanePoint, "");
  cla.AddArgument("refpoint", A_M_SPAC, &argReferencePoint, "");
  cla.AddArgument("applytransl", A_SPAC, &applytransl, "");
  if (!cla.Parse())
    return false; // bad arguments or wrong number of arguments ...

  // process found arguments:
  if (argInplanePoint.size() == 2)
  {
    for (unsigned int i = 0; i < argInplanePoint.size(); ++i)
      inplanepoint[i] = argInplanePoint[i];
    inplanepoint[2] = 0; // in zero-plane
  }
  else
    return false; // wrong number of arguments
  if (argReferencePoint.size() == 3)
  {
    for (unsigned int i = 0; i < argReferencePoint.size(); ++i)
      referencepoint[i] = argReferencePoint[i];
  }
  else
    return false; // wrong number of arguments
  
  // compute transformation:
  // WCS ... world coordinate system
  // F ... focal point in WCS
  // R ... reference point in WCS
  // R' ... reference point's projection in image-plane in WCS
  // X(x) ... x-coordinate of point X
  double transx, transy;
  InternalPointType origin = m_MovingImage->GetOrigin();
  typename BaseProjectionFilterType::InputPointType focalpoint =
    m_ProjectionFilter->GetFocalPoint();
  
  // -> add volume offset to reference point --> WCS!
  for (int z = 0; z < origin.PointDimension; z++)
    referencepoint[z] = referencepoint[z] + origin[z];
  // reference-point transformation in 'physical space' (WCS):
  referencepoint = m_Transformation->TransformPoint(referencepoint);
  
  // translation in x-direction: (F(x) - R'(x)) * R(z) / F(z) + R'(x) - R(x)
  transx = (focalpoint[0] - inplanepoint[0]) *
           referencepoint[2] / focalpoint[2] +
           inplanepoint[0] - referencepoint[0];
  // translation in y-direction: (F(y) - R'(y)) * R(z) / F(z) + R'(y) - R(y)
  transy = (focalpoint[1] - inplanepoint[1]) *
           referencepoint[2] / focalpoint[2] +
           inplanepoint[1] - referencepoint[1];

  if (applytransl)
  {
    // apply transformation as translation of transformation-component:
    typename BaseTransformType::OutputVectorType trans;

    trans = m_Transformation->GetTranslation();
    trans[0] += transx;
    trans[1] += transy;
    m_Transformation->SetTranslation(trans);

    m_Registration->SetInitialTransformParameters(
        m_Transformation->GetParameters()); // update initial transformation!
  }
  else
  {
    // apply transformation as origin of the volume:
    typename BaseTransformType::InputPointType cor;

    origin[0] += transx;
    origin[1] += transy;
    m_MovingImage->SetOrigin(origin);

    cor = m_Transformation->GetCenter();
    cor[0] += transx;
    cor[1] += transy;
    m_Transformation->SetCenter(cor); // update center of rotation!
  }

  return true;
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::StartRegistration() throw (
  ExceptionObject &)
{
  if (!m_Transformation || !m_Optimizer || !m_Registration ||
      !m_Interpolator || !m_ProjectionFilter || !m_Metric ||
      !m_FixedImage || !m_MovingImage)
  {
    itkExceptionMacro(<< "At least one essential registration component " <<
      "has not been connected, initialized or configured!");
    return;
  }

  // also offer the possibility to invoke registration by directly calling
  // StartRegistration() - therefore, check whether Update() or 
  // StartRegistration() has been called; Update() also calls 
  // StartRegistration(), but m_Updating is then expected to equal TRUE

  if (!this->m_Updating)
  {
    this->Update();
  }
  else
  {
    // store initial transformation:
    m_InitialTransformationParamters = 
      m_Registration->GetInitialTransformParameters();
    
    // Connect interpolator and internal perspective projection image filter
    // now:
    m_Interpolator->SetPerspectiveProjectionImageFilter(m_ProjectionFilter);
  
    // this->InvokeEvent(StartEvent()); // signal registration start (done by
    // Superclass)
  
    try
    {
      m_Registration->Update(); // trigger the pipeline
    }
    catch (ExceptionObject &e)
    {
      itkExceptionMacro(<< std::endl << 
        "An error occured during internal registration: " << e);
      this->InvokeEvent(AbortEvent()); // signal registration abort
      return;
    }
  
    // this->InvokeEvent(EndEvent()); // signal registration end (done by
    // Superclass)
  }
}

template <typename TInternalPixelType, typename TScalarType>
void 
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::StopRegistration()
{
  if (m_Registration)
    m_Registration->StopRegistration();
}

template <typename TInternalPixelType, typename TScalarType>
void 
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GenerateData()
{
  try
  {
    StartRegistration();
  }
  catch (ExceptionObject &e)
  {
    itkExceptionMacro(<< "ERROR DURING REGISTRATION: " << e);
  }
}

template <typename TInternalPixelType, typename TScalarType>
std::vector<unsigned long>
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::AddGenericFrameworkObserver(Command *observer) throw (ExceptionObject &)
{
  std::vector<unsigned long> cmdTags;

  if (!observer)
  {
    itkExceptionMacro(<< "Observer is NULL!");
    return cmdTags;
  }

  // start event
  cmdTags.push_back(this->AddObserver(StartEvent(), observer));
  // abort event
  cmdTags.push_back(this->AddObserver(AbortEvent(), observer));
  // end event
  cmdTags.push_back(this->AddObserver(EndEvent(), observer));

  return cmdTags;
}

template <typename TInternalPixelType, typename TScalarType>
unsigned long
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::AddMultiResolutionObserver(Command *observer) throw (ExceptionObject &)
{
  if (!m_Registration || !observer)
  {
    itkExceptionMacro(<< "Registration or observer is NULL!");
    return 0;
  }
  
  return m_Registration->AddObserver(IterationEvent(), observer);
}

template <typename TInternalPixelType, typename TScalarType>
unsigned long
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::AddOptimizationObserver(Command *observer) throw (ExceptionObject &)
{
  if (!m_Registration || !m_Optimizer || !observer)
  {
    itkExceptionMacro(<< "Registration or optimizer or observer is NULL!");
    return 0;
  }
  
  m_OptimizerObservers.push_back(observer);
  
  return m_Optimizer->AddObserver(IterationEvent(), observer);
}

template <typename TInternalPixelType, typename TScalarType>
void 
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::ReAddOptimizationObservers()
{
	if (!m_Registration || !m_Optimizer)
  {
    itkExceptionMacro(<< "Registration or optimizer is NULL!");
    return;
  }
  
  for (unsigned int i = 0; i < m_OptimizerObservers.size(); i++)
  	m_Optimizer->AddObserver(IterationEvent(), m_OptimizerObservers[i]);
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::RescaleFixedImageInRegion(InternalPixelType minIntensity, 
  InternalPixelType maxIntensity)
{
  typedef itk::ImageRegionConstIterator<InternalImageType> InputIteratorType;
  typedef itk::ImageRegionIterator<InternalImageType> OutputIteratorType;
  typedef itk::CastImageFilter<InternalImageType, InternalImageType>
    CastFilterType;
  typedef typename CastFilterType::Pointer CastFilterPointer;
  
  InternalSpacingType fixedSpacing;
  InternalRegionType fixedRegion;
  InternalSpacingType movingSpacing;  
  
  // make a copy of the fixed image:
  CastFilterPointer castFilter = CastFilterType::New();
  
  castFilter->SetInput(m_FixedImage);
  m_RescaledFixedImage = castFilter->GetOutput();
  
  try
  {
    castFilter->Update();
  }
  catch (ExceptionObject &e)
  {
    std::cerr << "ERROR copying image: " << e << std::endl;
  }
  m_RescaledFixedImage->DisconnectPipeline();
  
  InputIteratorType iit(m_RescaledFixedImage, m_FixedImageRegion);
  InternalPixelType realMin = NumericTraits<InternalPixelType>::max(); 
  InternalPixelType realMax = NumericTraits<InternalPixelType>::Zero;
  InternalPixelType v;
  double scale, offset;
  
  // find current min / max intensity values in fixed image region:
  for (iit.GoToBegin(); !iit.IsAtEnd(); ++iit)
  {
    v = iit.Get();
    
    if (v > realMax)
      realMax = v;
    if (v < realMin)
      realMin = v;
  }
  
  // modify the copied fixed image:
  OutputIteratorType oit(m_RescaledFixedImage, m_FixedImageRegion);
  
  // scale and shift:
  if (realMax != realMin)
    scale = (maxIntensity - minIntensity) / (realMax - realMin);
  else if (realMax != 0)
    scale = (maxIntensity - minIntensity) / realMax;
  else
    scale = 0.;
  offset = minIntensity - realMin * scale;

  for (oit.GoToBegin(); !oit.IsAtEnd(); ++oit)
  {
    v = oit.Get();        
    
    // apply rescaling:
    v = v * scale + offset;                
    
    oit.Set(v);
  }   
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GetCurrentImageProperties(InternalSpacingType &fixedSpacing,
  InternalRegionType &fixedRegion, InternalSpacingType &movingSpacing)
{
  if (!m_Registration || !m_Registration->GetFixedImagePyramid() ||
      !m_Registration->GetMovingImagePyramid())
    return;

  int level = m_Registration->GetCurrentLevel();
  ScheduleType fsched = m_Registration->GetFixedImagePyramid()->GetSchedule();
  InternalSpacingType fspac = m_FixedImage->GetSpacing();

  for (unsigned int d = 0; d < InternalDimension; d++)
    fixedSpacing[d] = static_cast<typename InternalSpacingType::ValueType>(
      static_cast<float>(fsched[level][d]) *
      static_cast<float>(fspac[d]));

  ScheduleType msched = m_Registration->GetMovingImagePyramid()->GetSchedule();
  InternalSpacingType mspac = m_MovingImage->GetSpacing();

  for (unsigned int d = 0; d < InternalDimension; d++)
    movingSpacing[d] = static_cast<typename InternalSpacingType::ValueType>(
        static_cast<float>(msched[level][d]) *
        static_cast<float>(mspac[d]));

  // compute current fixed image region:
  typename InternalRegionType::SizeType fsize = m_FixedImageRegion.GetSize();
  typename InternalRegionType::IndexType fidx = m_FixedImageRegion.GetIndex();
  typename InternalRegionType::SizeType size;
  typename InternalRegionType::IndexType idx;
  for (unsigned int d = 0; d < InternalDimension; d++)
  {
    const float sf = static_cast<float>(fsched[level][d]);

    size[d] = static_cast<typename InternalRegionType::SizeType::SizeValueType>(
      vcl_floor(static_cast<float>(fsize[d]) / sf));
    if(size[d] < 1)
      size[d] = static_cast<typename InternalRegionType::SizeType::
        SizeValueType>(1);

    idx[d] = static_cast<typename InternalRegionType::IndexType::
      IndexValueType>(vcl_ceil(static_cast<float>(fidx[d]) / sf));
  }
  fixedRegion.SetSize(size);
  fixedRegion.SetIndex(idx);
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::AdjustProjectorOutputImageProperties()
{
  if (!m_Registration || !m_Registration->GetFixedImagePyramid() ||
      !m_Registration->GetMovingImagePyramid() ||
      !m_ProjectionFilter || !m_Interpolator)
    return;

  // origin of projection-filter remains the same (in world coordinate system)

  // adjust X-ray-region to resolution level:

  typename BaseProjectionFilterType::OutputSizeType outsize;
  typename InternalRegionType::SizeType origsize = m_XrayImageRegion.GetSize();
  InternalSpacingType outspac;

  ScheduleType mps = m_Registration->GetMovingImagePyramid()->GetSchedule();
  unsigned int level = m_Registration->GetCurrentLevel();

  for (unsigned int d = 0; d < InternalDimension; ++d)
  {
    outsize[d] =
      static_cast<typename InternalRegionType::SizeType::SizeValueType>(
      static_cast<float>(origsize[d]) /
      static_cast<float>(mps[level][d]));
    if(outsize[d] < 1)
      outsize[d] = static_cast<typename InternalRegionType::SizeType::
        SizeValueType>(1);

    outspac[d] = static_cast<float>(m_XrayImageSpacing[d] * mps[level][d]);
  }

  m_ProjectionFilter->SetOutputSize(outsize);
  m_ProjectionFilter->SetOutputSpacing(outspac);
}

template <typename TInternalPixelType, typename TScalarType>
typename MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType,
  TScalarType>::InternalImageConstPointer
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GetCurrentProjectionImage()
{
  if (!m_MovingImage || !m_ProjectionFilter || !m_Transformation ||
      !m_Interpolator)
   return NULL;

  try
  {
    return m_Interpolator->GetProjectedImage();
  }
  catch (ExceptionObject &e)
  {
    std::cerr << "Could not retrieve current moving image with current " <<
      "configuration: " << e << std::endl;
    return NULL;
  }
}

template <typename TInternalPixelType, typename TScalarType>
typename MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType,
  TScalarType>::HistImagePointer
MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, TScalarType>
::GetCurrentMetricHistogram(bool update)
{

  /* typedefs for concrete casting */
  typedef NormalizedMutualInformationHistogramImageToImageMetric<
    InternalImageType, InternalImageType> NmiMetricType;
  typedef typename NmiMetricType::Pointer NmiMetricPointer;

  NmiMetricPointer nmiDummy = NmiMetricType::New();

  if (m_Metric && (m_Metric->GetNameOfClass() == nmiDummy->GetNameOfClass()))
  {
    typedef typename NmiMetricType::HistogramType HistogramType;
    typedef HistogramToEntropyImageFilter<HistogramType>
      HistogramToEntropyImageFilterType;
    typedef Image<double, 2> HistFilterImageType;
    //typedef CastImageFilter<HistFilterImageType, HistImageType>
    //  CastFilterType;
    typedef RescaleIntensityImageFilter<HistFilterImageType, HistImageType>
      RescaleFilterType;

    NmiMetricPointer metric = SmartPointer<NmiMetricType>(
      static_cast<NmiMetricType *>(m_Metric.GetPointer()));
    typename HistogramToEntropyImageFilterType::Pointer histEntropyFilter =
      HistogramToEntropyImageFilterType::New();
    //CastFilterType::Pointer castFilter = CastFilterType::New();
    RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();    

    if (update)
      metric->GetValue(metric->GetTransform()->GetParameters());

    histEntropyFilter->SetInput(metric->GetHistogram());
    //castFilter->SetInput(histEntropyFilter->GetOutput());
    rescaleFilter->SetInput(histEntropyFilter->GetOutput());
    rescaleFilter->SetOutputMinimum(0.);
    rescaleFilter->SetOutputMaximum(255.);

    //HistImageType::Pointer outImg = castFilter->GetOutput();
    //castFilter->Update();
    HistImageType::Pointer outImg = rescaleFilter->GetOutput();
    rescaleFilter->Update();

    // ensure that spacing is set to 1/1:
    HistImageType::SpacingType space;
    space[0] = 1.0;
    space[1] = 1.0;
    outImg->SetSpacing(space);
    outImg->DisconnectPipeline();

    return outImg;
  }
  else
    return NULL;
}


}


#endif /* ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_TXX_ */
