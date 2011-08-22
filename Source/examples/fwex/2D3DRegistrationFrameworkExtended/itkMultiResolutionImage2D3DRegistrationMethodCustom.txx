

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_TXX_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_TXX_


#include <itkMultiResolutionImage2D3DRegistrationMethodCustom.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkMath.h>

// support the configure methods:
#include <itkMattesMutualInformationImageToImageMetricComplete.h>


namespace itk
{


template <typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, TScalarType>
::MultiResolutionImage2D3DRegistrationMethodCustom()
  : Superclass()
{
  ;
}

template <typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, TScalarType>
::~MultiResolutionImage2D3DRegistrationMethodCustom()
{
  ;
}

template <typename TInternalPixelType, typename TScalarType>
void
MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, TScalarType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << " -- Custom 2D/3D Registration Framework -- " << std::endl;

}

template <typename TInternalPixelType, typename TScalarType>
bool
MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, TScalarType>
::ConfigureMetric(const std::string config)
{
  if (!this->m_Metric)
    return false;

  // type casting:
  typedef MattesMutualInformationImageToImageMetricComplete
    <InternalImageType, InternalImageType> MmiMetricType;
  typedef typename MmiMetricType::Pointer MmiMetricPointer;
  
  MmiMetricPointer mmiMetric = MmiMetricType::New(); // dummy

  // obviously it is a custom-configuration (MMI):
  if (this->m_Metric->GetNameOfClass() == mmiMetric->GetNameOfClass()) 
  {     
    itksys::CommandLineArguments cla;
    int argMmiBins = -1;
    itksys_stl::vector<double> argMmiLb;
    itksys_stl::vector<double> argMmiUb;
    double argNoOfSamples = -1; // default: all
    unsigned long noOfSamples = 1;

    // parse configure string:
    this->CallArgsInitialize(config, &cla);
    cla.AddArgument("histbins", this->A_SPAC, &argMmiBins, "");
    cla.AddArgument("lowerbounds", this->A_M_SPAC, &argMmiLb, "");
    cla.AddArgument("upperbounds", this->A_M_SPAC, &argMmiUb, "");
    cla.AddArgument("nosamples", this->A_SPAC, &argNoOfSamples, "");
    if (!cla.Parse())
      return false; // bad arguments or wrong number of arguments ...

    // process found arguments:
    if (argMmiBins <= 0)
      return false; // need number of bins

    // set-up metric:
    mmiMetric = SmartPointer<MmiMetricType>(static_cast<MmiMetricType *>(
      this->m_Metric.GetPointer()));            
    
    if (argNoOfSamples < 0)
    {
      // NOTE: there is a difference between UseAllPixels and 
      // NumberOfSpatialSamples>=image-pixel-count>!
      mmiMetric->SetUseAllPixels(true);
    }
    else
    { 
      if (argNoOfSamples > 1.0) // discrete pixel count
      {
        noOfSamples = itk::Math::Round<unsigned long, double>(argNoOfSamples);        
      }
      else // relative pixel count (0.0 .. 1.0)
      {
        typename Superclass::InternalSpacingType fixedSpacing;
        typename Superclass::InternalRegionType fixedRegion;
        typename Superclass::InternalSpacingType movingSpacing;
        
        // extract current fixed image region information:
        this->GetCurrentImageProperties(fixedSpacing, fixedRegion, 
          movingSpacing);
        
        unsigned long total = fixedRegion.GetSize()[0] *
          fixedRegion.GetSize()[1];
          
        if (total > 0)
        {                    
          noOfSamples = itk::Math::Round<unsigned long, double>(
            static_cast<double>(total) * argNoOfSamples);
        }
        else
          return false; // region is invalid
      }
      
      mmiMetric->SetUseAllPixels(false);
      mmiMetric->SetNumberOfSpatialSamples(noOfSamples);
    }
    
    mmiMetric->SetNumberOfHistogramBins(argMmiBins);
    
    mmiMetric->SetLowerBound(argMmiLb); // if size() == 0/1/2 is internally checked
    mmiMetric->SetUpperBound(argMmiUb);

    return true;
  }
  else
    return Superclass::ConfigureMetric(config); // delegate
}

template <typename TInternalPixelType, typename TScalarType>
typename MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, 
  TScalarType>::HistImagePointer
MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType, TScalarType>
::GetCurrentMetricHistogram(bool update)
{
  // type casting:
  typedef MattesMutualInformationImageToImageMetricComplete<
    InternalImageType, InternalImageType> MmiMetricType;
  typedef typename MmiMetricType::Pointer MmiMetricPointer;
  
  MmiMetricPointer mmiMetric = MmiMetricType::New(); // dummy

  // Mattes Mutual Information:
  if (this->m_Metric && 
      this->m_Metric->GetNameOfClass() == mmiMetric->GetNameOfClass()) 
  { 
    typedef typename MmiMetricType::JointPDFType JointPDFType;
    
    mmiMetric = SmartPointer<MmiMetricType>(static_cast<MmiMetricType *>(
      this->m_Metric.GetPointer()));

    if (update) // update demanded
    { 
      mmiMetric->GetValue(mmiMetric->GetTransform()->GetParameters());
    }

    typename JointPDFType::Pointer pdf = mmiMetric->GetJointPDF();

    if (pdf) // -> rescale to 0..255
    { 
      typedef RescaleIntensityImageFilter<JointPDFType, JointPDFType>
	RescaleFilterType;

      typename RescaleFilterType::Pointer rescaler = RescaleFilterType::New();

      rescaler->SetInput(pdf);
      rescaler->SetOutputMinimum(0.);
      rescaler->SetOutputMaximum(255.);

      typename JointPDFType::Pointer outImg = rescaler->GetOutput();

      rescaler->Update();

      outImg->DisconnectPipeline();

      return outImg;
    }
    else
      return NULL;
  }
  else // delegate
    return Superclass::GetCurrentMetricHistogram(update);
}


}


#endif /* ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_TXX_ */
