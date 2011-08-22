

#ifndef ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_H_
#define ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_H_


#include <itkCommand.h>

#include <itkMultiResolutionRegistrationCommand.h>

#include <fstream>


namespace itk
{


/** \class MultiResolutionRegistrationCommandCustom
 *  \brief This is a customized observer with minimal implementation for 
 *  observing multi-resolution level changes of the extended internal 
 *  registration.
 *
 * <p>This support class provides a simple implementation of a multi-resolution
 * registration observer.
 *
 * @see itk::MultiResolutionRegistrationCommand
 * 
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.0
 */
template <typename TFramework>
class MultiResolutionRegistrationCommandCustom
  : public MultiResolutionRegistrationCommand<TFramework>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionRegistrationCommandCustom Self;
  typedef MultiResolutionRegistrationCommand<TFramework> Superclass;
  typedef SmartPointer<Self> Pointer;

  itkNewMacro(Self);

  /** basic types **/
  typedef TFramework FrameworkType;
  typedef typename FrameworkType::Pointer FrameworkPointer;
  typedef typename FrameworkType::RegistrationType RegistrationType;
  typedef typename RegistrationType::Pointer RegistrationPointer;
  typedef typename FrameworkType::InternalImageType InternalImageType;
  typedef typename InternalImageType::Pointer InternalImagePointer;

  /** metric types **/
  typedef typename FrameworkType::BaseMetricType BaseMetricType;
  typedef typename BaseMetricType::Pointer BaseMetricPointer;

protected:
  /** Constructor of the command. **/
  MultiResolutionRegistrationCommandCustom();

  /** Destructor of the command. **/
  virtual ~MultiResolutionRegistrationCommandCustom();
  
  /**
   * Implement re-configuration of registration components (new metric: MMI)
   * level-changes (>0).
   * @see Superclass::ReConfigureRegistrationComponents()
   */  
  virtual void ReConfigureRegistrationComponents();
  
};


}


#include "itkMultiResolutionRegistrationCommandCustom.txx"

#endif /* ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_H_ */
