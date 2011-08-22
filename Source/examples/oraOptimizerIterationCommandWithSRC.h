

#ifndef ORAOPTIMIZERITERATIONCOMMANDWITHSRC_H_
#define ORAOPTIMIZERITERATIONCOMMANDWITHSRC_H_


#include <itkOptimizerIterationCommand.h>


namespace ora /** open radART **/
{


/** \class OptimizerIterationCommandWithSRC
 *  \brief This is a base class with minimal implementation for observing
 *  iteration events of the internal optimizer.
 *
 * <p>This support class provides a simple implementation of an optimizer
 * iteration observer. It could be used as base-class for customized
 * optimizer observers in combination with
 * <code>ora::MultiResolutionImage2D3DRegistrationMethodWithSRC</code>.</p>
 *
 * @see ora::MultiResolutionImage2D3DRegistrationMethodWithSRC
 * @see itk::OptimizerIterationCommand
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.1
 */
template <typename TFramework>
class OptimizerIterationCommandWithSRC
  : public itk::OptimizerIterationCommand<TFramework>
{
public:
  /** Standard class typedefs. */
  typedef OptimizerIterationCommandWithSRC Self;
  typedef itk::OptimizerIterationCommand<TFramework> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

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


protected:
  /** Optimizer iteration observer constructor. **/
  OptimizerIterationCommandWithSRC();
  /** Optimizer iteration observer destructor. **/
  virtual ~OptimizerIterationCommandWithSRC();

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
  
  /**
   * Write sample distribution at current optimization iteration (SRC).
   * @param currIt current optimization iteration
   * @param costValueChanged flag indicating that cost-value has changed since
   * last iteration
   */
  virtual void WritePDFImage(int currIt, bool costValueChanged);

};


}


#include "oraOptimizerIterationCommandWithSRC.txx"


#endif /* ORAOPTIMIZERITERATIONCOMMANDWITHSRC_H_ */
