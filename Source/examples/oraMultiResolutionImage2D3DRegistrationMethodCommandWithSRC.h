

#ifndef ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_H_
#define ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_H_


#include <itkMultiResolutionImage2D3DRegistrationMethodCommand.h>

namespace ora /** open radART **/
{

/** \class MultiResolutionImage2D3DRegistrationMethodCommandWithSRC
 * \brief A base class with minimal implementation for observing the basic
 * events of MultiResolutionImage2D3DRegistrationMethodWithSRC.
 *
 * <p>This support class provides a simple implementation of an
 * MultiResolutionImage2D3DRegistrationMethodWithSRC observer. It could be
 * used as base-class for customized registration observers.</p>
 *
 * In contrast to its superclass, this class implements additional support
 * for stochastic rank correlation metric, regular step gradient descent
 * optimization and writing of "sample distribution" files (instead of PDFs if
 * SRC is the current metric).
 *
 * @see itk::MultiResolutionImage2D3DRegistrationMethodCommand
 * @see ora::MultiResolutionImage2D3DRegistrationMethodWithSRC
 *
 * @author phil <phil.steininger e_T gmail.com>
 * @version 1.1
 */
template <typename TFramework>
class MultiResolutionImage2D3DRegistrationMethodCommandWithSRC
  : public itk::MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionImage2D3DRegistrationMethodCommandWithSRC Self;
  typedef itk::MultiResolutionImage2D3DRegistrationMethodCommand<TFramework> Superclass;
  typedef itk::SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;

protected:
  /** Constructor of the command. **/
  MultiResolutionImage2D3DRegistrationMethodCommandWithSRC();

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

  /**
   * Write out final sample distribution (SRC). Different interpretation for
   * SRC than for other metrics!
   * @param fw pointer to 2D/3D registration framework
   * @param finalPars final transformation parameters
   */
  virtual void WriteFinalPDFImage(FrameworkPointer fw,
    typename FrameworkType::BaseOptimizerType::ParametersType finalPars);

};


}


#include "oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.txx"


#endif /* ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMANDWITHSRC_H_ */
