

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_H_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_H_


#include <itkCommand.h>
#include <fstream>


namespace itk
{


/** \class MultiResolutionImage2D3DRegistrationMethodCommand
 *  \brief This is a base class with minimal implementation for observing
 *  the basic events of MultiResolutionImage2D3DRegistrationMethod.
 *
 * <p>This support class provides a simple implementation of an
 * MultiResolutionImage2D3DRegistrationMethod observer. It could be used as 
 * base-class for customized registration observers.</p>
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.2
 */
template <typename TFramework>
class MultiResolutionImage2D3DRegistrationMethodCommand
  : public Command
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionImage2D3DRegistrationMethodCommand Self;
  typedef Command Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;

  /** Setter/getter for verbose-output-mode. **/
  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  /**
   * Setter/getter for log file name (file specification).
   * If Logging is set to ON and file name is specified the final registration
   * parameters (last level) are logged!
   * NOTE: a "%d" can be applied in the file name where the final level
   * number is inserted.
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
   */
  itkSetMacro(Logging, bool);
  itkGetMacro(Logging, bool);

  /**
   * Setter/getter for image-output mode (true=ON, false=OFF).
   * If image-output is set to ON and the 2 image file names are set the
   * images are automatically generated and written!
   */
  itkSetMacro(ImageOutput, bool);
  itkGetMacro(ImageOutput, bool);

  /**
   * Setter/getter for image file names in the following order:
   *  final moving image <br>
   *  final PDF image. <br>
   *
   * If image-output is set to ON and the 2 image file names are set the
   * images are automatically generated and written!
   *
   * NOTE: you can also add an empty string ("") in order to leave the according
   * image out, e.g. fileNames[0]="", fileNames[1]="/path/finalPDF.mhd" to
   * force final PDF output, but leave out final moving image output!
   */
  virtual void SetImageFileNames(std::vector<std::string> fileNames)
  { 
    if (fileNames.size() >= 2)
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

protected:
  /** verbose-output-mode **/
  bool m_Verbose;
  /** log file name **/
  std::string m_LogFileName;
  /** logging mode **/
  bool m_Logging;
  /** log file itself **/
  std::ofstream m_LogFile;
  /** image-output-mode **/
  bool m_ImageOutput;
  /** image file names (final moving image, final PDF image) **/
  std::vector<std::string> m_ImageFileNames;
  /** time stamp of registration start **/
  double m_StartTimeStamp;

  /** Constructor of the command. **/
  MultiResolutionImage2D3DRegistrationMethodCommand();

  /** Specialized Execute()-method for start events. **/
  virtual void ExecuteStartEvent(Object *object);

  /** Specialized Execute()-method for abort events. **/
  virtual void ExecuteAbortEvent(Object *object);

  /** Specialized Execute()-method for end events. **/
  virtual void ExecuteEndEvent(Object *object);

  /**
   * Write out final moving image.
   * @param fw pointer to 2D/3D registration framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalMovingImage(FrameworkPointer fw,
    typename FrameworkType::BaseOptimizerType::ParametersType finalPars);
  /**
   * Write out final PDF image.
   * @param fw pointer to 2D/3D registration framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalPDFImage(FrameworkPointer fw,
    typename FrameworkType::BaseOptimizerType::ParametersType finalPars);

  /**
   * Verbose and log the final registration result.
   * @param fw pointer to 2D/3D registration framework
   * @param numberOfIterations number of optimization iterations
   * @param finalPars final transformation parameters
   * @param bestValue best metric value
   */
  virtual void VerboseAndLogFinalRegistrationResult(FrameworkPointer fw,
    unsigned int numberOfIterations,
    typename FrameworkType::BaseOptimizerType::ParametersType finalPars,
    double bestValue);

  /**
   * Get current registration results (optimizer/metric-dependent).
   * @param registration multi-resolution registration framework pointer
   * @param finalPars returned final transformation parameters
   * @param bestValue returned best metric value
   * @param numberOfIterations returned number of total iterations
   */
  virtual bool GetCurrentResults(typename FrameworkType::Pointer fw,
    typename FrameworkType::BaseOptimizerType::ParametersType &finalPars,
    double &bestValue, unsigned int &numberOfIterations);

};


}


#include "itkMultiResolutionImage2D3DRegistrationMethodCommand.txx"

#endif /* ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_H_ */
