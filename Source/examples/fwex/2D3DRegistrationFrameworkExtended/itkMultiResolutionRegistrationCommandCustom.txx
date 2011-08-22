

#ifndef ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_TXX_
#define ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_TXX_

#include <itkMultiResolutionRegistrationCommandCustom.h>

// support for concrete casting:
#include <itkMattesMutualInformationImageToImageMetricComplete.h>


namespace itk
{


template <typename TFramework>
MultiResolutionRegistrationCommandCustom<TFramework>
::MultiResolutionRegistrationCommandCustom()
	: MultiResolutionRegistrationCommand<TFramework>()
{
	;
}

template <typename TFramework>
MultiResolutionRegistrationCommandCustom<TFramework>
::~MultiResolutionRegistrationCommandCustom()
{
  ;
}

template <typename TFramework>
void
MultiResolutionRegistrationCommandCustom<TFramework>
::ReConfigureRegistrationComponents()
{  
  RegistrationPointer registration = this->m_Framework->GetRegistration();
  unsigned int level = registration->GetCurrentLevel(); 
  if (level < 1)
    return;

  // metric
  if (this->m_MetricSels.size() >= level && 
      this->m_MetricConfigs.size() >= this->m_MetricSels.size()) // need config
  {
  	if (this->m_MetricSels[level - 1] == "MMI")
  	{
	    if (this->m_Verbose)
	      std::cout << std::endl << "Resetting metric component for level " <<
	        (level + 1) << ": " << this->m_MetricSels[level - 1];    
	        
      typedef MattesMutualInformationImageToImageMetricComplete<
        InternalImageType, InternalImageType> MmiMetricType;
      typename MmiMetricType::Pointer mmi = MmiMetricType::New();
      mmi->SetFixedImage(registration->GetFixedImage());
      mmi->SetMovingImage(registration->GetMovingImage());
      mmi->SetTransform(registration->GetTransform());
      mmi->SetInterpolator(registration->GetInterpolator());
      typename InternalImageType::RegionType freg = registration->GetMetric()->
        GetFixedImageRegion();
      mmi->SetFixedImageRegion(freg);      
      registration->SetMetric(mmi);
      this->m_Framework->SetMetric(registration->GetMetric());
      
	    if (this->m_Verbose)
	    	std::cout << " --> SUCCESS" << std::endl;
	  }
  }
 	// metric configuration is handled via base-class!

	// superclass implementation:
  this->Superclass::ReConfigureRegistrationComponents();
}


}


#endif /* ITKMULTIRESOLUTIONREGISTRATIONCOMMANDCUSTOM_TXX_ */
