

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_H_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_H_


#include <itkProcessObject.h>
#include <itkImage.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkMultiResolutionPyramidImageFilter.h>
#include <itkMatrixOffsetTransformBase.h>
#include <itkImageToImageMetric.h>
#include <itkSingleValuedNonLinearOptimizer.h>
#include <itkRealTimeClock.h>
#include <itkDataObjectDecorator.h>
#include <itksys/CommandLineArguments.hxx>

#include <itkMultiResolutionImage2D3DRegistrationMethodCommand.h>
#include <itkMultiResolutionRegistrationCommand.h>
#include <itkOptimizerIterationCommand.h>
#include <itkExhaustiveOptimizerComplete.hxx>
#include <itkPerspectiveProjectionInterpolateImageFunction.h>
#include <itkPerspectiveProjectionImageFilter.h>

#include <vector>


namespace itk
{


/** \class MultiResolutionImage2D3DRegistrationMethod
 *  \brief This class encapsulates the 2D/3D registration framework for X-rays.
 *
 * <p>This class simply encapsulates a typical 2D/3D registration framework
 * compatible with splatting and raycast interpolators.</p>
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 2.0
 */
template <typename TInternalPixelType = float, typename TScalarType = double>
class MultiResolutionImage2D3DRegistrationMethod
  : public ProcessObject
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionImage2D3DRegistrationMethod Self;
  typedef ProcessObject Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /* basic typedefs */
  /** image type applied internally **/
  typedef TInternalPixelType InternalPixelType;
  const static unsigned int InternalDimension = 3;
  typedef Image<InternalPixelType, InternalDimension> InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef typename InternalImageType::ConstPointer InternalImageConstPointer;
  typedef typename InternalImageType::RegionType InternalRegionType;
  typedef typename InternalImageType::PointType InternalPointType;
  typedef typename InternalImageType::SpacingType InternalSpacingType;
  /** histogram image type (metrics) **/
  typedef float HistPixelType;
  const static unsigned int HistDimension = 2;
  typedef Image<HistPixelType, HistDimension> HistImageType;
  typedef typename HistImageType::Pointer HistImagePointer;

  /** registration and components **/

  /** registration **/
  typedef MultiResolutionImageRegistrationMethod<InternalImageType,
    InternalImageType> RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;
  typedef MultiResolutionPyramidImageFilter<InternalImageType,
    InternalImageType> ImagePyramidType;
  typedef typename ImagePyramidType::Pointer ImagePyramidPointer;
  typedef typename ImagePyramidType::ScheduleType ScheduleType;

  /** transformation **/
  typedef TScalarType ScalarType;
  typedef MatrixOffsetTransformBase<ScalarType, 3, 3> BaseTransformType;
  typedef typename BaseTransformType::Pointer BaseTransformPointer;
  typedef typename BaseTransformType::ParametersType
    BaseTransformParametersType;

  /** interpolator **/
  typedef double CoordinateRepresentationType;
  typedef Collector::SumCollector<InternalPixelType> CollectorType;
  typedef PerspectiveProjectionInterpolateImageFunction<InternalImageType,
    CoordinateRepresentationType, CollectorType> BaseInterpolatorType;
  typedef typename BaseInterpolatorType::Pointer BaseInterpolatorPointer;
  typedef typename BaseInterpolatorType::PerspectiveProjectionImageFilterType
    BaseProjectionFilterType;
  typedef typename BaseProjectionFilterType::Pointer
    BaseProjectionFilterPointer;

  /** metric **/
  typedef ImageToImageMetric<InternalImageType, InternalImageType>
    BaseMetricType;
  typedef typename BaseMetricType::Pointer BaseMetricPointer;

  /** optimizer **/
  typedef SingleValuedNonLinearOptimizer BaseOptimizerType;
  typedef BaseOptimizerType::Pointer BaseOptimizerPointer;

  /** observers **/
  typedef MultiResolutionImage2D3DRegistrationMethodCommand<Self>
    GenericFrameworkObserverType;
  typedef typename GenericFrameworkObserverType::Pointer
    GenericFrameworkObserverPointer;
  typedef MultiResolutionRegistrationCommand<Self>
    MultiResolutionObserverType;
  typedef typename MultiResolutionObserverType::Pointer
    MultiResolutionObserverPointer;
  typedef OptimizerIterationCommand<Self> OptimizationObserverType;
  typedef typename OptimizationObserverType::Pointer
    OptimizationObserverPointer;

  /** clock tools **/
  typedef RealTimeClock::Pointer RealTimeClockPointer;
  typedef RealTimeClock::TimeStampType RealTimeStampType;

  /**
   * typedef for the output: using Decorator pattern for enabling
   * the Transform to be passed in the data pipeline
   */
  typedef DataObjectDecorator<BaseTransformType> BaseTransformOutputType;
  typedef typename BaseTransformOutputType::Pointer BaseTransformOutputPointer;
  typedef typename BaseTransformOutputType::ConstPointer
    BaseTransformOutputConstPointer;


  /** Run-time type information (and related methods). **/
  itkTypeMacro(MultiResolutionImage2D3DRegistrationMethod, ProcessObject);

  /** Method for creation through the object factory. **/
  itkNewMacro(Self);

  /** @return the transform resulting from the registration process */
  virtual const BaseTransformOutputType *GetOutput() const;

  /**
   * Make a DataObject of the correct type to be used as the specified
   * output.
   */
  virtual DataObjectPointer MakeOutput(unsigned int idx);

  /**
   * Method to return the latest modified time of this object or
   * any of its members
   */
  virtual unsigned long GetMTime() const;

  /** Set the fixed image with internal type. **/
  virtual void SetFixedImage(InternalImagePointer fixedImage)
  {
    m_FixedImage = fixedImage;
    this->Modified();
  }
  /** Get the fixed image with internal type. **/
  InternalImagePointer GetFixedImage()
  {
    if (m_FixedRescalingActivated)
      return m_RescaledFixedImage;
    else
      return m_FixedImage;
  }
  /** Set the fixed image region. **/
  virtual void SetFixedImageRegion(InternalRegionType fixedImageRegion)
  {
    m_FixedImageRegion = fixedImageRegion;
    m_FixedImageRegion.SetSize(2, 1);
    this->Modified();
  }
  /** Get the fixed image region. **/
  InternalRegionType GetFixedImageRegion()
  {
    return m_FixedImageRegion;
  }
  /** Set the moving image with internal type. **/
  virtual void SetMovingImage(InternalImagePointer movingImage)
  {
    m_MovingImage = movingImage;
    this->Modified();
  }
  /** Get the moving image with internal type. **/
  InternalImagePointer GetMovingImage()
  {
    return m_MovingImage;
  }
  /**
   * Set the X-ray image region (position and size of X-ray image plane).
   * NOTE: the 3rd dimension is ignored and is internally set to a static
   * index=0 and size=1!
   */
  virtual void SetXrayImageRegion(InternalRegionType xrayImageRegion)
  {
    m_XrayImageRegion = xrayImageRegion;
    m_XrayImageRegion.SetIndex(2, 0); // totally planar
    m_XrayImageRegion.SetSize(2, 1);
    this->Modified();
  }
  /**
   * Get the X-ray image region (position and size of X-ray image plane).
   * NOTE: the 3rd dimension is ignored and is internally set to a static
   * index=0 and size=1!
   */
  InternalRegionType GetXrayImageRegion()
  {
    return m_XrayImageRegion;
  }
  /**
   * Set the X-ray image spacing (spacing of X-ray image plane pixel).
   * NOTE: the 3rd dimension is ignored and is internally set to a static
   * spacing=1.0!
   */
  virtual void SetXrayImageSpacing(InternalSpacingType xrayImageSpacing)
  {
    m_XrayImageSpacing = xrayImageSpacing;
    m_XrayImageSpacing[2] = 1.0; // totally planar
  }
  /**
   * Get the X-ray image spacing (spacing of X-ray image plane pixel).
   * NOTE: the 3rd dimension is ignored and is internally set to a static
   * spacing=1.0!
   */
  InternalSpacingType GetXrayImageSpacing()
  {
    return m_XrayImageSpacing;
  }

  /** Initialize the X-ray registration framework internally.
   * Create and connect the involved components.
   * <b>NOTE:</b> the metric/optimizer/transform/interpolator types
   * should be fixed before.
   * @return true if successful
   **/
  virtual bool Initialize();

  /**
   * Deinitialize the X-ray registration framework internally. Disconnect
   * the internal components.
   */
  virtual void Deinitialize();

  /** Get the internal ITK registration framework. **/
  virtual RegistrationPointer GetRegistration()
  {
    return m_Registration;
  }
  /** Get the internal registration transformation component. **/
  virtual BaseTransformPointer GetTransformation()
  {
    return m_Transformation;
  }
  /** Set the internal registration transformation component. **/
  virtual void SetTransformation(BaseTransformPointer transform)
  {
    if (this->m_Transformation != transform)
    {
      this->m_Transformation = transform;
      this->Modified();
    }
  }
  /** Get the internal registration interpolation (projection) component. **/
  virtual BaseInterpolatorPointer GetInterpolator()
  {
    return m_Interpolator;
  }
  /** Get the internal projector component. **/
  virtual BaseProjectionFilterPointer GetProjectionFilter()
  {
    return m_ProjectionFilter;
  }
  /** Set the internal projector component. **/
  virtual void SetProjectionFilter(BaseProjectionFilterPointer filter)
  {
    if (this->m_ProjectionFilter != filter)
    {
      this->m_ProjectionFilter = filter;
      this->Modified();
    }
  }
  /** Get the internal registration metric component. **/
  virtual BaseMetricPointer GetMetric()
  {
    return m_Metric;
  }
  /** Set the internal registration metric component. **/
  virtual void SetMetric(BaseMetricPointer metric)
  {
    if (this->m_Metric != metric)
    {
      this->m_Metric = metric;
      this->Modified();
    }
  }
  /** Get the internal registration optimizer component. **/
  virtual BaseOptimizerPointer GetOptimizer()
  {
    return m_Optimizer;
  }
  /** Set the internal registration optimizer component. **/
  virtual void SetOptimizer(BaseOptimizerPointer optimizer)
  {
    if (this->m_Optimizer != optimizer)
    {
      this->m_Optimizer = optimizer;
      this->Modified();
    }
  }

  /** Get the optimization direction (min:TRUE/max:FALSE); depends on metric. **/
  virtual bool GetOptimizationMinimization()
  {
    return m_OptimizationMinimization;
  }

  /** Get activate/deactivate fixed image rescaling in fixed image region **/
  virtual bool GetFixedRescalingActivated()
  {
    return m_FixedRescalingActivated;
  }
  /** Set activate/deactivate fixed image rescaling in fixed image region **/
  virtual void SetFixedRescalingActivated(bool activate)
  {
    m_FixedRescalingActivated = activate;
    this->Modified();
  }
  /** Get minimum intensity for fixed imgage rescaling **/
  virtual InternalPixelType GetFixedRescalingMin()
  {
    return m_FixedRescalingMin;
  }
  /** Set minimum intensity for fixed imgage rescaling **/
  virtual void SetFixedRescalingMin(InternalPixelType intensity)
  {
    m_FixedRescalingMin = intensity;
    this->Modified();
  }
  /** Get maximum intensity for fixed imgage rescaling **/
  virtual InternalPixelType GetFixedRescalingMax()
  {
    return m_FixedRescalingMax;
  }
  /** Set maximum intensity for fixed imgage rescaling **/
  virtual void SetFixedRescalingMax(InternalPixelType intensity)
  {
    m_FixedRescalingMax = intensity;
    this->Modified();
  }

  /**
   * Get current image properties for current multi-resolution level according
   * to set image pyramids.
   * @param fixedSpacing returned current spacing of fixed image
   * @param fixedRegion returned current region of fixed image
   * @param movingSpacing returned current spacing of moving image
   */
  virtual void GetCurrentImageProperties(InternalSpacingType &fixedSpacing,
    InternalRegionType &fixedRegion, InternalSpacingType &movingSpacing);

  /**
   * Adjust projector's output image properties to current multi-resolution
   * level settings. <b>It is important to call that method when resolution
   * level changes occur!</b>
   */
  virtual void AdjustProjectorOutputImageProperties();

  /**
   * Get current projection image (depends on current projection and
   * transformation properties.
   */
  virtual InternalImageConstPointer GetCurrentProjectionImage();

  /**
   * Get metric's current histogram (only possible if metric is histogram-based,
   * e.g. NMI, MMI). The intensities of the histogram are rescaled to [0.,255.].
   * @param update if TRUE then metric->GetValue() is called to ensure that the
   * histogram corresponds to current transformation settings; otherwise the
   * metric is assumed to be up-to-date!
   */
  virtual HistImagePointer GetCurrentMetricHistogram(bool update);

  /**
   * Converts argument-string (space-separated) into typical argc/argv-style
   * and then initializes the specified command line parser.
   * @param s argument-string containing space-separated arguments
   * @param cla ITK command line parser
   */
  void CallArgsInitialize(std::string s, itksys::CommandLineArguments *cla);

//  /**
//   *
//   * NOTE: (for all the ConfigureXXX()-methods)
//   * - you should apply all possible arguments in the configuration-string in
//   *   order to achieve a transparent configuration of the components!
//   * - the order of arguments in configuration-string does not matter
//   * - arguments and values are separate by spaces, e.g. "-arg value"
//   * - "<type>"-symbols describe placeholders for values of a specified type
//   * - "..."-symbols describe a list of (semantically clear) values (multi-
//   *   arguments)
//   * - multi-arguments are space separated, e.g. "-arg val1 val2 val3"
//   * - applying unknown arguments results in failure of the ConfigureXXX()-
//   *   method!
//   *
//   */

  /**
   * Configure the set images by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * (MOVING IMAGE, FIXED IMAGE)
   * "volumeorigin <double> <double> <double> fixedrescaling <double> <double>"
   *
   * where
   *
   * (MOVING IMAGE)
   * <b>volumeorigin</b> (doubles) defines the 3D volume offset (origin),
   * (FIXED IMAGE)
   * <b>fixedrescaling</b> (doubles) activates the intensity rescaling of the
   * fixed image within its fixed image region (first value is desired minimum
   * intensity, second value is desired maximum intensity; second value must
   * be greater than first value)
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureImages(std::string config);

  /**
   * Configure the internal image pyramids by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * Level 0 has original spacing and all further levels are consecutively
   * scaled by factor 2.
   * @param config string expression of format:
   *
   * "levels <int>"
   *
   * where
   *
   * <b>levels</b> (int) represents the number of levels of the image
   * pyramid
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureImagePyramids(const std::string config);

  /**
   * Configure the internal transformation by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * (EULER3D)
   * "cor <double> <double> <double> rot <double> <double> <double>
   * transl <double> <double> <double>"
   *
   * (AFFINE)
   * "cor <double> <double> <double> rot <double> <double> <double>
   * transl <double> <double> <double> scale <double> <double> <double>"
   *
   * where
   *
   * NOTE:
   * general order: rotation(, scaling,) translation
   * rotation: rotation x, rotation y, rotation z around center of rotation
   * scaling: based on center of rotation!
   *
   * <b>cor</b> (doubles) represents the 3D center of rotation of the volume
   * (relative to the volume's origin),
   * <b>rot</b> (doubles) is the 3D rotation of the volume (in degrees) and
   * <b>transl</b> (doubles) is the 3D translation of the volume,
   * <b>scale</b> (doubles) is the scaling in each direction of the volume.
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureTransformation(const std::string config);

  /**
   * Configure the internal optimizer by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * (EXHAUSTIVE)
   * "min <bool> scales <double> <double> ... steplen <double>
   * stepcount <int> <int> ..."
   *
   * (EVOLUTIONARY)
   * "min <bool> scales <double> <double> ... iterations <int> nvgseed <int>
   * radius <double> gfact <double> sfact <double> epsilon <double>"
   *
   * where
   *
   * <b>min</b> (bool) specifies the optimization criterion (0=maximize,
   * 1=minimize),
   * <b>scales</b> (doubles) specifies the weightings of the
   * transformation dimensions,
   * <b>steplen</b> (double) specifies the
   * basic step length (which is weighted by scales in each dimension) and
   * <b>stepcount</b> (integers) specifies the number of steps for each dimension,
   * <b>iterations</b> (integer) specifies the maximum number of iteration to be
   * performed (unless radius falls below epsilon),
   * <b>nvgseed</b> (integer) specifies the random seed for the normal variate
   * generator,
   * <b>radius</b> (double) specifies the initial search radius,
   * <b>gfact</b> (double) specifies the growth factor for search radius,
   * <b>sfact</b> (double) specifies the shrink factor for search radius,
   * <b>epsilon</b> (double) specifies the minimum search radius.
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureOptimizer(const std::string config);

  /**
   * Configure the internal metric by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * (NMI)
   * "histbins <int> <int> lowerbounds <double> <double>
   * upperbounds <double> <double>"
   *
   * where
   *
   * <b>histbins</b> (integers) specifies the size of the histogram
   * (moving, fixed),
   * <b>lowerbounds</b> (doubles) specifies the lower intensity bounds of
   * the histogram (fixed, moving),
   * <b>upperbounds</b> (doubles) specifies the upper intensity bounds of the
   * histogram (fixed, moving).
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureMetric(const std::string config);

  /**
   * Configure the internal interpolator (projection component) by applying a
   * configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * (SPLATTING)
   * "jobs <int> threshold <double> transferfunc <int> <double> <double> ...
   * transferinterp <int> focalpoint <double> <double> <double>
   * wobbmode <int> wobbsigma <double>
   * wobbcount <int> splatpixassign <int> wobbdeterm <bool>
   * rescale <double> <double>"
   *
   * (RAYCASTING)
   * "jobs <int> threshold <double> transferfunc <int> <double> <double> ...
   * transferinterp <int> focalpoint <double> <double> <double>
   * rescale <double> <double> raystepsize <double> raystepinterp <int>"
   *
   * where
   *
   * <b>jobs</b> (integer) specifies the number of threads (if jobs is <=0 then
   * the filter automatically detects the number of CPUs),
   * <b>threshold</b> (double) is the primary intensity threshold for the voxels
   * of interest,
   * <b>transferfunc</b> (integer, doubles) is specifies the transfer-function
   * (1stelement is the number of pairs, next elements are the input-output-
   * pairs),
   * <b>transferinterp</b> (integer) specifies the transfer-function
   * interpolation method (0:NONE 1:STEPWISE 2:LINEAR),
   * <b>focalpoint</b> (doubles) specifies the 3D-coordinate of the focalpoint
   * in world coordinate system in mm,
   * <b>rescale</b> (doubles) specifies an optional intermediate rescale-
   * intensity-filter which is applied to the projected image right after its
   * generation (thefirst argument specifies the minimum and the second the
   * maximum intensity),
   * <b>wobbmode</b> (integer) selects the focal offset computation (wobbling)
   * method (0:NONE 1:PRE_CONTINUOUS 2:PRE_RANDOM 3:RANDOM),
   * <b>wobbsigma</b> (doubles) specifies the focal offset standard deviation,
   * <b>wobbcount</b> (integer) specifies the number of pre-computed focal
   * offset values,
   * <b>splatpixassign</b> (integer) selects the projected pixel assignment
   * method to its neighbours (0:NEAREST 1:NEIGHBOURS),
   * <b>wobbdeterm</b> (bool 0/1) determines whether focal offset is computed
   * deterministically (reproducable) or not,
   * <b>raystepsize</b> (double) sets step size between interpolation points,
   * <b>raystepinterp</b> (integer) determines the interpolation method (
   * intensity at ray-points) (0:NearestNeighbour 1:Linear 2:BSpline(Order 3)
   * 3:WindowedSinc(Welch, Window 5))
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureInterpolator(const std::string config);

  /**
   * Computes the necessary x/y-transformation of the volume in order to
   * align the in-plane-projection of the reference point (relative to
   * volume) with the desired in-plane-position (x/y)! The resulting x/y-
   * translation is either applied as a x/y-translation within the
   * transformation component or a volume offset transformation of the volume.
   * <b>It is important that focal point, volume offset and moving image have
   * already been set when this method is invoked!</b>
   * (This method should simplify the process of initializing the
   * transformation, sth. like a rule-of-thumb intialization.)
   * @param config string expression of format:
   *
   * "inplane <double> <double> refpoint <double> <double> <double>
   * applytransl <bool>"
   *
   * where
   *
   * <b>inplane</b> (doubles) specifies the desired x/y-position of the
   * projection of the reference point in mm in world coordinate system (z-
   * coordinate of image-plane is defined as 0),
   * <b>refpoint</b> (doubles) specifies the 3D-coordinate of a reference point
   * relative to the volume's origin in mm (not in world coordinate system),
   * <b>applytransl</b> (bool 0/1) 0:apply resulting x/y-transformation as
   * translation of the transformation-component, 1:apply resulting x/y-
   * transformation as translation of the volume origin.
   *
   * @result true if configuration string has successfully been applied and
   * the translation has successfully been computed
   */
  virtual bool ComputeVolumeTranslationByReferencePoint(
    const std::string config);

  /**
   * @deprecated
   * Does the registration with the connected and configured components.
   * @throws itk::ExceptionObject if any any errors occur.
   */
  virtual void StartRegistration() throw (ExceptionObject &);

  /**
   * Stops the registration.
   */
  virtual void StopRegistration();

  /**
   * Add a generic observer (manages several commands) for monitoring the
   * MultiResolutionImage2D3DRegistrationMethod.
   * @param observer the generic monitoring observer object
   * @return a vector of command-tags that can be used for deleting or
   * retrieving the generic commands
   * @throws itk::ExceptionObject if any errors occur.
   */
  virtual std::vector<unsigned long> AddGenericFrameworkObserver(
    Command *observer) throw (ExceptionObject &);

  /**
   * Add an observer for monitoring the internal multi-resolution level
   * changes.
   * NOTE: the internal registration component must be initialized before
   * adding the observer!
   * @param observer the monitoring observer object
   * @return the command-tag that can be used for deleting or retrieving the
   * command from the internal components
   * @throws itk::ExceptionObject if any errors occur.
   */
  virtual unsigned long AddMultiResolutionObserver(
    Command *observer) throw (ExceptionObject &);

  /**
   * Add an observer for monitoring the optimization iterations.
   * NOTE: the internal registration and optimizer components must be
   * initialized before adding the observer!
   * @param observer the monitoring observer object
   * @return the command-tag that can be used for deleting or retrieving the
   * command from the internal components
   * @throws itk::ExceptionObject if any errors occur.
   */
  virtual unsigned long AddOptimizationObserver(
    Command *observer) throw (ExceptionObject &);

  /** Get reference to internal clock object **/
  RealTimeClockPointer GetClock()
  {
    return m_Clock;
  }

  /** Get last internal time stamp. **/
  RealTimeStampType GetLastTimeStamp()
  {
    return m_LastTimeStamp;
  }
  /** Set last internal time stamp. **/
  void SetLastTimeStamp(RealTimeStampType timeStamp)
  {
    m_LastTimeStamp = timeStamp;
  }

  /** intial transformation parameters **/
  BaseTransformParametersType GetInitialTransformationParamters()
  {
    return m_InitialTransformationParamters;
  }

  /**
   * Add the optimization observers again (important when exhanging
   * components.
   **/
  virtual void ReAddOptimizationObservers();

protected:
  /** shorter constants for shorter code **/
  static const itksys::CommandLineArguments::ArgumentTypeEnum A_SPAC =
      itksys::CommandLineArguments::SPACE_ARGUMENT;
  static const itksys::CommandLineArguments::ArgumentTypeEnum A_M_SPAC =
      itksys::CommandLineArguments::MULTI_ARGUMENT;

  /** fixed image **/
  InternalImagePointer m_FixedImage;
  /** rescaled fixed image **/
  InternalImagePointer m_RescaledFixedImage;
  /** fixed image region **/
  InternalRegionType m_FixedImageRegion;
  /** moving image **/
  InternalImagePointer m_MovingImage;
  /** X-ray image region (position and size of X-ray image plane) **/
  InternalRegionType m_XrayImageRegion;
  /** X-ray image spacing (spacing of X-ray image plane pixels) **/
  InternalSpacingType m_XrayImageSpacing;

  /** ITK registration framework (multi-resolutional) **/
  RegistrationPointer m_Registration;
  /** registration transformation **/
  BaseTransformPointer m_Transformation;
  /** registration interpolate image function **/
  BaseInterpolatorPointer m_Interpolator;
  /** internal projection filter (core of interpolator) **/
  BaseProjectionFilterPointer m_ProjectionFilter;
  /** registration intensity-based metric **/
  BaseMetricPointer m_Metric;
  /** registration intensity-based metric **/
  BaseOptimizerPointer m_Optimizer;

  /** specifies optimization direction (depends on selected metric) **/
  bool m_OptimizationMinimization;

  /** activate/deactivate fixed image rescaling in fixed image region **/
  bool m_FixedRescalingActivated;
  /** minimum intensity for fixed imgage rescaling **/
  InternalPixelType m_FixedRescalingMin;
  /** maximum intensity for fixed imgage rescaling **/
  InternalPixelType m_FixedRescalingMax;

  /** intial transformation parameters **/
  BaseTransformParametersType m_InitialTransformationParamters;

  /** just a helper **/
  RealTimeClockPointer m_Clock;
  /** last time stamp - can be used as central memory **/
  RealTimeStampType m_LastTimeStamp;

  /** stores added optimizer observers **/
  std::vector<Command *> m_OptimizerObservers;

  /** Hidden constructor. **/
  MultiResolutionImage2D3DRegistrationMethod();

  /** Discarded destructor. **/
  virtual ~MultiResolutionImage2D3DRegistrationMethod();

  /** Print class information. **/
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

  /**
   * Rescale the fixed image in the fixed image region to a specified interval.
   * @param minIntensity desired minimum intensity of the rescaled image
   * @param maxIntensity desired maximum intensity of the rescaled image
   **/
  virtual void RescaleFixedImageInRegion(InternalPixelType minIntensity,
    InternalPixelType maxIntensity);

  /**
   * Method invoked by the pipeline in order to trigger the computation of
   * the registration.
   */
  virtual void GenerateData();


private:
  // purposely not implemented
  MultiResolutionImage2D3DRegistrationMethod(const Self&);
  // purposely not implemented
  void operator=(const Self&);
};


}


#include "itkMultiResolutionImage2D3DRegistrationMethod.txx"


#endif /*ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHOD_H_*/
