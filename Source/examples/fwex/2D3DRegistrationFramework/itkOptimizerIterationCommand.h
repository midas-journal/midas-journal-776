

#ifndef ITKOPTIMIZERITERATIONCOMMAND_H_
#define ITKOPTIMIZERITERATIONCOMMAND_H_


#include <itkCommand.h>
#include <itkRealTimeClock.h>

#include <fstream>


namespace itk
{


/** \class OptimizerIterationCommand
 *  \brief This is a base class with minimal implementation for observing
 *  iteration events of the internal optimizer.
 *
 * <p>This support class provides a simple implementation of an optimizer
 * iteration observer. It could be used as base-class for customized
 * optimizer observers in combination with
 * <code>itk::MultiResolutionImage2D3DRegistrationMethod</code>.</p>
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.2
 */
template <typename TFramework>
class OptimizerIterationCommand
  : public Command
{
public:
  /** Standard class typedefs. */
  typedef OptimizerIterationCommand Self;
  typedef Command Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro(Self);


  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;

  /** optimizer types **/
  typedef typename FrameworkType::BaseOptimizerType BaseOptimizerType;
  typedef typename BaseOptimizerType::Pointer BaseOptimizerPointer;
  typedef typename BaseOptimizerType::ConstPointer BaseOptimizerConstPointer;

  /** transform types **/
  typedef typename FrameworkType::BaseTransformType BaseTransformType;
  typedef typename BaseTransformType::Pointer BaseTransformPointer;

  /** metric types **/
  typedef typename FrameworkType::BaseMetricType BaseMetricType;
  typedef typename BaseMetricType::Pointer BaseMetricPointer;

  /** Setter/getter for verbose-output-mode. **/
  itkSetMacro(Verbose, bool);
  itkGetMacro(Verbose, bool);

  /**
   * Setter/getter for log file name (file specification). Log file name can
   * contain a "%d" where the current level is automatically inserted!
   * If Logging is set to ON and file name specified the transform parameters
   * and metric values are logged!
   */
  virtual void SetLogFileName(std::string logFileName)
  {
    m_LogFileName = logFileName;
  }
  virtual std::string GetLogFileName()
  { return m_LogFileName;
  }
  /**
   * Setter/getter for logging mode (true=ON, false=OFF).
   * If Logging is set to ON and file name specified the transform parameters
   * and metric values are logged!
   */
  itkSetMacro(Logging, bool);
  itkGetMacro(Logging, bool);

  /**
   * Setter/getter for current parameter logging mode (true=ON, false=OFF).
   * If this mode is switched on current parameters (current transformation
   * parameters and resulting metric value) are append to the ordinary
   * logging parameters. <br>
   * <b>NOTE:</b> this mode slows down registration process because the
   * metric value unfortunately must be called twice!
   */
  itkSetMacro(LogCurrentParametersAlso, bool);
  itkGetMacro(LogCurrentParametersAlso, bool);

  /**
   * Setter/getter for log split value (log file is split into log files with
   * n lines; -1 for no splitting).
   */
  virtual void SetLogSplitValue(int n)
  {
    m_LogSplitValue = n;
  }
  virtual int GetLogSplitValue()
  {
    return m_LogSplitValue;
  }


  /** Setter/getter for image-output-mode. **/
  itkSetMacro(ImageOutput, bool);
  itkGetMacro(ImageOutput, bool);

  /**
   * Setter/getter for automatic-image-output-mode. That means that whenever
   * the COST VALUE CHANGES an image is written out (provided that
   * ImageBaseFileName is set). The ImageModulo attribute is also considered -
   * should be set <=0 to avoid confusions. ImageOutput must be set to TRUE
   * to enable ImageAutoOutput! PLEASE NOTE that this mode makes no sense
   * for exhaustive optimization because the cost value is likely to change
   * each iteration ...
   */
  itkSetMacro(ImageAutoOutput, bool);
  itkGetMacro(ImageAutoOutput, bool);

  /**
   * Setter/getter for image base file name. An image is generated each
   * n-th iteration (ImageModulo) and written to disc if image-output-mode
   * is active and ImageBaseFileName is set to a valid file name. In the
   * file name specification an optional %d can be applied which is a
   * placeholder for the current multi-resolution level. Before the last "."
   * the current iteration number is inserted.
   */
  virtual void SetImageBaseFileName(std::string fileName)
  {
    m_ImageBaseFileName = fileName;
  }
  virtual std::string GetImageBaseFileName()
  {
    return m_ImageBaseFileName;
  }

  /**
   * Setter/getter for image modulo value (an image is generated each n-th
   * iteration).
   */
  itkSetMacro(ImageModulo, unsigned int);
  itkGetMacro(ImageModulo, unsigned int);


  /** Setter/getter for PDF-output-mode (histogram). **/
  itkSetMacro(PDFOutput, bool);
  itkGetMacro(PDFOutput, bool);

  /**
   * Setter/getter for automatic-PDF-output-mode. That means that whenever
   * the COST VALUE CHANGES a PDF is written out (provided that
   * PDFBaseFileName is set). The PDFModulo attribute is also considered -
   * should be set <=0 to avoid confusions. PDFOutput must be set to TRUE
   * to enable PDFAutoOutput! PLEASE NOTE that this mode makes no sense
   * for exhaustive optimization because the cost value is likely to change
   * each iteration ...
   */
  itkSetMacro(PDFAutoOutput, bool);
  itkGetMacro(PDFAutoOutput, bool);

  /**
   * Setter/getter for PDF base file name. A PDF-image is generated each
   * n-th iteration (PDFModulo) and written to disc if PDF-output-mode
   * is active and PDFBaseFileName is set to a valid file name. In the
   * file name specification an optional %d can be applied which is a
   * placeholder for the current multi-resolution level. Before the last "."
   * the current iteration number is inserted.
   */
  virtual void SetPDFBaseFileName(std::string fileName)
  {
    m_PDFBaseFileName = fileName;
  }
  virtual std::string SetPDFBaseFileName()
  {
    return m_PDFBaseFileName;
  }

  /**
   * Setter/getter for PDF modulo value (a PDF-image is generated each n-th
   * iteration).
   */
  itkSetMacro(PDFModulo, unsigned int);
  itkGetMacro(PDFModulo, unsigned int);

  /**
   * Setter/getter for framework reference (needed for image-output-mode).
   */
  itkSetObjectMacro(Framework, FrameworkType);
  itkGetObjectMacro(Framework, FrameworkType);


  /** purposely not implemented. **/
  void Execute(Object *caller, const EventObject &event);

  /**
   * Main-method for observer.
   * @param object reference to the registration
   * @param event specifies the event-type
   */
  virtual void Execute(const Object *object, const EventObject &event);  

protected:
  /** verbose-output-mode **/
  bool m_Verbose;
  /** timing variables **/
  double m_StartTime;
  /** internal clock for time measurements **/
  itk::RealTimeClock::Pointer m_Clock;
  /** log file name **/
  std::string m_LogFileName;
  /** logging mode **/
  bool m_Logging;
  /** additional logging parameters **/
  bool m_LogCurrentParametersAlso;
  /** log file itself **/
  std::ofstream m_LogFile;
  /** splits the log file after N iterations **/
  int m_LogSplitValue;
  /** current log split ID **/
  int m_CurrLogSplitId;
  /** image-output-mode **/
  bool m_ImageOutput;
  /** automatic-image-output-mode **/
  bool m_ImageAutoOutput;
  /** image base filename (%d is optional level placeholder) **/
  std::string m_ImageBaseFileName;
  /** an image is generated each n-th iteration **/
  unsigned int m_ImageModulo;
  /** PDF-output-mode **/
  bool m_PDFOutput;
  /** automatic-PDF-output-mode **/
  bool m_PDFAutoOutput;
  /** PDF base filename (%d is optional level placeholder) **/
  std::string m_PDFBaseFileName;
  /** a PDF-image is generated each n-th iteration **/
  unsigned int m_PDFModulo;
  /** reference to the framework **/
  FrameworkPointer m_Framework;
  /** store cost value of previous iteration **/
  typename FrameworkType::BaseMetricType::MeasureType m_LastCostValue;

  /** Optimizer iteration observer constructor. **/
  OptimizerIterationCommand();

  /** Optimizer iteration observer destructor. **/
  virtual ~OptimizerIterationCommand();

  /**
   * Verbose registration progress information.
   * @param currIt current optimization iteration
   * @param maxIt maximum optimization iteration
   */
  virtual void VerboseProgress(int currIt, int maxIt);

  /**
   * Write out moving image at current optimization iteration.
   * @param currIt current optimization iteration
   * @param costValueChanged flag indicating that cost-value has changed since
   * last iteration
   */
  virtual void WriteMovingImage(int currIt, bool costValueChanged);
  /**
   * Write out PDF image at current optimization iteration.
   * @param currIt current optimization iteration
   * @param costValueChanged flag indicating that cost-value has changed since
   * last iteration
   */
  virtual void WritePDFImage(int currIt, bool costValueChanged);
  /**
   * Open/close/split central log file if necessary.
   * @param currIt current optimization iteration
   */
  virtual void OpenCloseSplitLogFile(int currIt);

  /**
   * Verbose and log the registration parameters (optimizer, transformation,
   * metric) at current iteration.
   * @param currIt current optimization iteration
   * @param thisIterationCostValue optimizer's cost value at this iteration
   * @param optimizerPos optimizer's position at this iteration
   */
  virtual void VerboseAndLogRegistrationParameters(int currIt,
    typename BaseMetricType::MeasureType thisIterationCostValue,
    typename BaseOptimizerType::ParametersType optimizerPos);

  /**
   * Get current iteration parameters depending on the type of optimizer.
   * @param optimizer reference to the optimizer
   * @param currIt returned current iteration number
   * @param maxIt returned maximum number of iterations
   * @param thisIterationCostValue returned cost-value of the current iteration
   * @param pos returned optimizer position of current iteration
   * @return true if the parameters could be retrieved
   */
  virtual bool GetCurrentIterationParameters(
    BaseOptimizerConstPointer optimizer,
    int &currIt, int &maxIt,
    typename BaseMetricType::MeasureType &thisIterationCostValue,
    typename BaseOptimizerType::ParametersType &pos);
    
  /**
   * Stops the optimizer. This is important for immediate stopping of the 
   * registration (multi-resolution framework's stop-method shows effect at
   * resolution-level-changes only).
   */
  virtual void StopOptimizer();
  
};


}


#include "itkOptimizerIterationCommand.txx"


#endif /* ITKOPTIMIZERITERATIONCOMMAND_H_ */
