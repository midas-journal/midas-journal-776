

#ifndef ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_H_
#define ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_H_


#include <itkCommand.h>

#include <fstream>


namespace itk
{


/** \class MultiResolutionRegistrationCommand
 *  \brief This is a base class with minimal implementation for observing
 *  multi-resolution level changes of the internal registration.
 *
 * <p>This support class provides a simple implementation of a multi-resolution
 * registration observer. It could be used as base-class for customized
 * registration observers in combination with
 * <code>itk::MultiResolutionImage2D3DRegistrationMethod</code>.</p>
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.4
 */
template <typename TFramework>
class MultiResolutionRegistrationCommand
  : public Command
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionRegistrationCommand Self;
  typedef Command Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;
  typedef typename FrameworkType::RegistrationType RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;
  typedef typename FrameworkType::InternalImageType InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef typename FrameworkType::InternalSpacingType InternalSpacingType;
  typedef typename FrameworkType::InternalRegionType InternalRegionType;

  /** optimizer types **/
  typedef typename FrameworkType::BaseOptimizerType BaseOptimizerType;
  typedef typename BaseOptimizerType::Pointer BaseOptimizerPointer;

  /** metric types **/
  typedef typename FrameworkType::BaseMetricType BaseMetricType;
  typedef typename BaseMetricType::Pointer BaseMetricPointer;

  /** Setter/getter for verbose-output-mode. **/
  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  /**
   * Setter/getter for log file name (file specification).
   * If Logging is set to ON and file name is specified the final registration
   * parameters of each (except the last) level are logged!
   */
  virtual void SetLogFileName(std::string logFileName)
  {
    m_LogFileName = logFileName;
  }
  virtual std::string GetLogFileName()
  {
    return m_LogFileName;
  }
  /**
   * Setter/getter for logging mode (true=ON, false=OFF).
   * If Logging is set to ON and file name specified the final parameters are
   * logged!
   * NOTE: a "%d" can be applied in the file name where the current level
   * number is inserted.
   */
  itkSetMacro(Logging, bool);
  itkGetMacro(Logging, bool);

  /**
   * Setter/getter for framework reference.
   */
  itkSetObjectMacro(Framework, FrameworkType);
  itkGetObjectMacro(Framework, FrameworkType);

  /**
   * Setter/getter for image-output mode (true=ON, false=OFF).
   * If image-output is set to ON and 3 (+2 optional) image file names are set
   * the images are automatically generated and written!
   */
  itkSetMacro(ImageOutput, bool);
  itkGetMacro(ImageOutput, bool);

  /**
   * Setter/getter for image file names in the following order:
   *  fixed image (should include "%d"-pattern) <br>
   *  input volume image (should include "%d"-pattern) <br>
   *  initial moving image (should include "%d"-pattern) <br>
   *  final moving image (should include "%d"-pattern) <br>
   *  final PDF image (should include "%d"-pattern). <br>
   *  cropped fixed image (should include "%d"-pattern). <br>
   *
   * If image-output is set to ON and the 5 image file names are
   * set the images are automatically generated and written!
   *
   * NOTE: you can also add an empty string ("") in order to leave the according
   * image out, e.g. fileNames[0]="", fileNames[1]="/path/inp-vol.mhd" ... to
   * force input volume output, but leave out fixed image output!
   */
  virtual void SetImageFileNames(std::vector<std::string> fileNames)
  {
    m_ImageFileNames = fileNames;
  }
  virtual std::vector<std::string> GetImageFileNames()
  {
    return m_ImageFileNames;
  }

  /**
   * Main-method for observer.
   * @param object reference to the registration
   * @param event specifies the event-type
   */
  virtual void Execute(Object *object, const EventObject &event);

  // purposely not implemented.
  void Execute(const Object *, const EventObject &);

  /**
   * Set metric selection strings for levels > 0. If not provided, the
   * metric component from previous level will remain. NOTE: for each 
   * metric selection we need a metric configuration string (but not vice versa)
   */
  virtual void SetMetricSels(const std::vector<std::string> &sels);
  /*
   * Get metric selection strings for levels > 0. If not provided, the
   * metric component from previous level will remain. NOTE: for each 
   * metric selection we need a metric configuration string (but not vice versa)
   */
  virtual const std::vector<std::string> &GetMetricSels() const;  
  /**
   * Set metric configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual void SetMetricConfigs(const std::vector<std::string> &configs);
  /**
   * Get metric configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual const std::vector<std::string> &GetMetricConfigs() const;

  /**
   * Set Optimizer selection strings for levels > 0. If not provided, the
   * Optimizer component from previous level will remain. NOTE: for each 
   * Optimizer selection we need a Optimizer configuration string (but not vice 
   * versa)
   */
  virtual void SetOptimizerSels(const std::vector<std::string> &sels);
  /*
   * Get Optimizer selection strings for levels > 0. If not provided, the
   * Optimizer component from previous level will remain. NOTE: for each 
   * Optimizer selection we need a Optimizer configuration string (but not vice 
   * versa)
   */
  virtual const std::vector<std::string> &GetOptimizerSels() const;  
  /**
   * Set Optimizer configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual void SetOptimizerConfigs(const std::vector<std::string> &configs);
  /**
   * Get Optimizer configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual const std::vector<std::string> &GetOptimizerConfigs() const;
  
  /**
   * Set Interpolator configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual void SetInterpolatorConfigs(const std::vector<std::string> &configs);
  /**
   * Get Interpolator configuration strings for levels > 0. If not provided, the
   * configuration from previous level will remain.
   */
  virtual const std::vector<std::string> &GetInterpolatorConfigs() const;

protected:
  /** verbose-output-mode **/
  bool m_Verbose;
  /** logging mode ON or OFF **/
  bool m_Logging;
  /** log file's name specification **/
  std::string m_LogFileName;
  /** log file itself **/
  std::ofstream m_LogFile;
  /** pointer to the MultiResolutionImage2D3DRegistrationMethod **/
  FrameworkPointer m_Framework;
  /** image output mode ON/OFF **/
  bool m_ImageOutput;
  /** file names of the output image files **/
  std::vector<std::string> m_ImageFileNames;
  /** metric selection strings for levels > 0 **/
  std::vector<std::string> m_MetricSels;
  /** metric configuration strings for levels > 0 **/
  std::vector<std::string> m_MetricConfigs;
  /** optimizer selection strings for levels > 0 **/
  std::vector<std::string> m_OptimizerSels;
  /** optimizer configuration strings for levels > 0 **/
  std::vector<std::string> m_OptimizerConfigs;
  /** interpolator configuration strings for levels > 0 **/
  std::vector<std::string> m_InterpolatorConfigs;

  /** Constructor of the command. **/
  MultiResolutionRegistrationCommand();

  /** Destructor of the command. **/
  virtual ~MultiResolutionRegistrationCommand();

  /**
   * Write out the fixed image (if demanded).
   * @param registration reference to the multi resolution framework
   */
  virtual void WriteFixedImage(RegistrationPointer registration);
  /**
   * Write out the input volume (if demanded).
   * @param registration reference to the multi resolution framework
   */
  virtual void WriteInputVolume(RegistrationPointer registration);
  /**
   * Write out the moving image (if demanded).
   * @param registration reference to the multi resolution framework
   */
  virtual void WriteMovingImage(RegistrationPointer registration);
  /**
   * Write out the final moving image (if demanded).
   * @param registration reference to the multi resolution framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalMovingImage(RegistrationPointer registration,
    typename BaseOptimizerType::ParametersType finalPars);
  /**
   * Write out the final PDF image (if demanded).
   * @param registration reference to the multi resolution framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalPDFImage(RegistrationPointer registration,
    typename BaseOptimizerType::ParametersType finalPars);

  /**
   * Simple output of the current spacing settings (moving, fixed ...)
   * due to the image pyramid settings.
   * @param registration reference to the multi resolution framework
   */
  virtual void VerboseSpacingProperties(RegistrationPointer registration);

  /**
   * Summarized output of registration result at current level.
   * @param registration reference to the multi resolution framework
   * @param numberOfIterations number of iterations at level
   * @param finalPars the final registration parameters
   * @param bestValue best value of metric (depends on min/max-mode)
   */
  virtual void VerboseAndLogRegistrationResultAtLevel(
    RegistrationPointer registration, unsigned int numberOfIterations,
    typename BaseOptimizerType::ParametersType finalPars, double bestValue);

  /**
   * Get current registration results (optimizer/metric-dependent).
   * @param registration multi-resolution registration framework pointer
   * @param finalPars returned final transformation parameters
   * @param bestValue returned best metric value
   * @param numberOfIterations returned number of total iterations
   */
  virtual bool GetCurrentResults(RegistrationPointer registration,
    typename BaseOptimizerType::ParametersType &finalPars,
    double &bestValue, unsigned int &numberOfIterations);

  /**
   * In some special cases (e.g. exhaustive optimizer) the initial position
   * for next the next level must be set manually.
   */
  virtual void TakeOverInitialPositionForNextLevel(
    RegistrationPointer registration);
  
  /**
   * Implement re-configuration of registration components (optimizer, metric,
   * interpolator) during level-changes (>0).
   */  
  virtual void ReConfigureRegistrationComponents();
  
};


}


#include "itkMultiResolutionRegistrationCommand.txx"

#endif /* ITKMULTIRESOLUTIONREGISTRATIONCOMMAND_H_ */
