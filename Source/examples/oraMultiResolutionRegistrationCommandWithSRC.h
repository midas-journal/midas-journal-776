

#ifndef ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_H_
#define ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_H_


#include <itkMultiResolutionRegistrationCommandCustom.h>


namespace ora /** open radART **/
{


/** \class MultiResolutionRegistrationCommandWithSRC
 *  \brief This is an extended base class with minimal implementation for
 *  observing multi-resolution level changes of the internal registration
 *  compatible with MultiResolutionImage2D3DRegistrationMethodCommandWithSRC.
 *
 * <p>This support class provides a simple implementation of a multi-resolution
 * registration observer. It could be used as base-class for customized
 * registration observers in combination with <code>
 * ora::MultiResolutionImage2D3DRegistrationMethodCommandWithSRC</code>.</p>
 *
 * NOTE: The call to SetImageFileNames() can now take one more file name:
 * fixed image mask file (should include "%d"-pattern). <br>
 * @see itk::MultiResolutionRegistrationCommand#SetImageFileNames()
 *
 * @see ora::MultiResolutionImage2D3DRegistrationMethodCommandWithSRC
 * @see itk::MultiResolutionRegistrationCommandCustom
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.1
 */
template <typename TFramework>
class MultiResolutionRegistrationCommandWithSRC
  : public itk::MultiResolutionRegistrationCommandCustom<TFramework>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionRegistrationCommandWithSRC Self;
  typedef itk::MultiResolutionRegistrationCommandCustom<TFramework> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;
  typedef typename FrameworkType::RegistrationType RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;
  typedef typename FrameworkType::InternalImageType InternalImageType;
  typedef typename FrameworkType::RankImageType RankImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;
  typedef typename FrameworkType::InternalSpacingType InternalSpacingType;
  typedef typename FrameworkType::InternalRegionType InternalRegionType;

  /** optimizer types **/
  typedef typename FrameworkType::BaseOptimizerType BaseOptimizerType;
  typedef typename BaseOptimizerType::Pointer BaseOptimizerPointer;

  /** metric types **/
  typedef typename FrameworkType::BaseMetricType BaseMetricType;
  typedef typename BaseMetricType::Pointer BaseMetricPointer;

  /** @see itk::MultiResolutionRegistrationCommand::Execute() **/
  virtual void Execute(itk::Object *object, const itk::EventObject &event);

protected:

  /** Constructor of the command. **/
  MultiResolutionRegistrationCommandWithSRC();

  /** Destructor of the command. **/
  virtual ~MultiResolutionRegistrationCommandWithSRC();

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
   * Implement re-configuration of registration components (optimizer, metric,
   * interpolator) during level-changes (>0).
   */  
  virtual void ReConfigureRegistrationComponents();
  
  /**
   * Write out the final sample distribution (SRC).
   * @param registration reference to the multi resolution framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalPDFImage(RegistrationPointer registration,
    typename BaseOptimizerType::ParametersType finalPars);

  /** Write fixed image mask for current level if defined. **/
  virtual void WriteFixedImageMask();

};


}


#include "oraMultiResolutionRegistrationCommandWithSRC.txx"


#endif /* ORAMULTIRESOLUTIONREGISTRATIONCOMMANDWITHSRC_H_ */
