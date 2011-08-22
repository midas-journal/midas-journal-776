
#ifndef ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_TXX_
#define ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_TXX_

//
#include "oraMultiResolutionImage2D3DRegistrationMethodWithSRC.h"

#include "oraStochasticRankCorrelationImageToImageMetric.h"

#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkImageMaskSpatialObject.h>

namespace ora /** open radART **/
{

template<typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
    TScalarType>::MultiResolutionImage2D3DRegistrationMethodWithSRC() :
  Superclass()
{

}

template<typename TInternalPixelType, typename TScalarType>
MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
    TScalarType>::~MultiResolutionImage2D3DRegistrationMethodWithSRC()
{

}

template<typename TInternalPixelType, typename TScalarType>
void MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
    TScalarType>::PrintSelf(std::ostream& os, itk::Indent indent) const
{
  this->Superclass::PrintSelf(os, indent);

  os << indent << " -- 2D/3D Registration Framework with SRC -- " << std::endl;
}

template<typename TInternalPixelType, typename TScalarType>
bool MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
    TScalarType>::ConfigureMetric(const std::string config)
{
  if (!this->m_Metric)
    return false;

  // type casting:
  typedef StochasticRankCorrelationImageToImageMetric<InternalImageType,
      InternalImageType, RankImageType> SRCMetricType;
  typedef typename SRCMetricType::Pointer SRCMetricPointer;

  SRCMetricPointer srcMetric = SRCMetricType::New(); // dummy

  bool retval;

  // obviously it is a custom-configuration (SRC):
  if (this->m_Metric->GetNameOfClass() == srcMetric->GetNameOfClass())
  {
    itksys::CommandLineArguments cla;
    std::vector<int> argSRCRandomSeeds;
    std::vector<double> argSRCDerivativeScales;
    std::vector<double> argSRCFixedHistogram;
    std::vector<double> argSRCMovingHistogram;
    double argSRCSampleCoverage = -1;
    bool argSRCHorn = false;
    bool argSRC0Ranks = true;

    // parse configure string:
    this->CallArgsInitialize(config, &cla);
    cla.AddArgument("rseeds", this->A_M_SPAC, &argSRCRandomSeeds, "");
    cla.AddArgument("dscales", this->A_M_SPAC, &argSRCDerivativeScales, "");
    cla.AddArgument("fhist", this->A_M_SPAC, &argSRCFixedHistogram, "");
    cla.AddArgument("mhist", this->A_M_SPAC, &argSRCMovingHistogram, "");
    cla.AddArgument("coverage", this->A_SPAC, &argSRCSampleCoverage, "");
    cla.AddArgument("horn", this->A_SPAC, &argSRCHorn, "");
    cla.AddArgument("zerorankscontribute", this->A_SPAC, &argSRC0Ranks, "");
    if (!cla.Parse())
      return false; // bad arguments or wrong number of arguments ...

    // cast!
    srcMetric = NULL;
    srcMetric = itk::SmartPointer<SRCMetricType>(static_cast<SRCMetricType *>(
      this->m_Metric.GetPointer()));

    typename SRCMetricType::SeedsType rseeds;
    rseeds.SetSize(3);
    rseeds.Fill(0);
    if (argSRCRandomSeeds.size() == 2)
    {
      rseeds[0] = argSRCRandomSeeds[0];
      rseeds[1] = argSRCRandomSeeds[1];
    } // else: default
    srcMetric->SetRandomSeeds(rseeds);

    typename SRCMetricType::ScalesType dscales;
    dscales.SetSize(this->m_Transformation->GetNumberOfParameters());
    dscales.Fill(1);
    if (argSRCDerivativeScales.size()
        == this->m_Transformation->GetNumberOfParameters())
    {
      for (unsigned int i = 0; i
          < this->m_Transformation->GetNumberOfParameters(); i++)
        dscales[i] = argSRCDerivativeScales[i];
    } // else: default
    srcMetric->SetDerivativeScales(dscales);

    // require fhist!
    if (argSRCFixedHistogram.size() == 4)
    {
      srcMetric->SetFixedHistogramMinIntensity(argSRCFixedHistogram[0]);
      srcMetric->SetFixedHistogramMaxIntensity(argSRCFixedHistogram[1]);
      srcMetric->SetFixedNumberOfHistogramBins(
          static_cast<int> (argSRCFixedHistogram[2]));
      srcMetric->SetFixedHistogramClipAtEnds(
          static_cast<bool> (argSRCFixedHistogram[3]));
    }
    else
      return false;

    // require mhist!
    if (argSRCMovingHistogram.size() == 4)
    {
      srcMetric->SetMovingHistogramMinIntensity(argSRCMovingHistogram[0]);
      srcMetric->SetMovingHistogramMaxIntensity(argSRCMovingHistogram[1]);
      srcMetric->SetMovingNumberOfHistogramBins(
          static_cast<int> (argSRCMovingHistogram[2]));
      srcMetric->SetMovingHistogramClipAtEnds(
          static_cast<bool> (argSRCMovingHistogram[3]));
    }
    else
      return false;

    if (argSRCSampleCoverage <= 0)
      return false;
    if (argSRCSampleCoverage > 100)
      argSRCSampleCoverage = 100;
    srcMetric->SetSampleCoverage(argSRCSampleCoverage);

    srcMetric->SetUseHornTiedRanksCorrection(argSRCHorn);
    srcMetric->SetMovingZeroRanksContributeToMeasure(argSRC0Ranks);

    srcMetric->SetNoOverlapReactionMode(1);
    srcMetric->SetNoOverlapMetricValue(1000);

    retval = true;
  }
  else
  {
    retval = Superclass::ConfigureMetric(config); // delegate
  }

  // check the new option which is applicable to ALL metrics (fixed mask):
  std::string::size_type p = config.find("circularfixedmask");
  if (this->m_FixedImage && p != std::string::npos)
  {
    double midp[2];
    double radius = 0;
    std::string s = config.substr(p + 18);
    int c = 0;
    while (c < 3)
    {
      p = s.find(" ");
      if (p != std::string::npos)
      {
        if (c < 2)
          midp[c] = atof(s.substr(0, p).c_str());
        else
          radius = atof(s.substr(0, p).c_str());
        s = s.substr(p + 1);
        c++;
      }
      else if (c == 2 && s.length() > 0) // last
      {
        radius = atof(s.c_str());
        c++;
      }
      else
      {
        break;
      }
    }
    if (c == 3 && radius > 0) // valid values -> generate mask
    {
      FixedMaskImagePointer fim = FixedMaskImageType::New();
      typename InternalImageType::RegionType freg =
          this->m_FixedImage->GetLargestPossibleRegion();
      fim->SetRegions(freg);
      fim->Allocate();
      fim->SetOrigin(this->m_FixedImage->GetOrigin());
      fim->SetSpacing(this->m_FixedImage->GetSpacing());
      fim->SetDirection(this->m_FixedImage->GetDirection());
      fim->FillBuffer(0);

      typedef itk::ImageRegionIterator<FixedMaskImageType> MaskIteratorType;
      MaskIteratorType fimit(fim, freg);
      double spac[2];
      spac[0] = this->m_FixedImage->GetSpacing()[0];
      spac[1] = this->m_FixedImage->GetSpacing()[1];
      radius = radius * radius; // sqr
      while (!fimit.IsAtEnd())
      {
        typename FixedMaskImageType::IndexType idx = fimit.GetIndex();
        double posX = (double)idx[0] * spac[0];
        double posY = (double)idx[1] * spac[1];

        double xx = posX - midp[0];
        double yy = posY - midp[1];
        if ((xx * xx + yy * yy) <= radius)
          fimit.Set(1); // set pixel inside circular mask
        ++fimit;
      }

      typedef itk::ImageMaskSpatialObject<InternalImageType::ImageDimension>
        MaskSpatialObjectType;
      typedef typename MaskSpatialObjectType::Pointer MaskSpatialObjectPointer;
      MaskSpatialObjectPointer fspatial = MaskSpatialObjectType::New();
      fspatial->SetImage(fim);
      fspatial->Update();

      // set the fixed image mask to the metric:
      this->m_Metric->SetFixedImageMask(fspatial);
    }
  }

  return retval;
}

template<typename TInternalPixelType, typename TScalarType>
bool MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
    TScalarType>::ConfigureOptimizer(const std::string config)
{
  typedef itk::RegularStepGradientDescentOptimizer RSGDOptimizerType;
  typedef typename RSGDOptimizerType::Pointer RSGDOptimizerPointer;

  RSGDOptimizerPointer rsgdOpt = RSGDOptimizerType::New(); // dummy

  // obviously it is a custom-configuration (RSGD):
  if (this->m_Optimizer->GetNameOfClass() == rsgdOpt->GetNameOfClass())
  {
    itksys::CommandLineArguments cla;
    int argRSGDMin = -1;
    std::vector<double> argRSGDScales;
    std::vector<double> argRSGDSteps;
    int argRSGDIterations = 100;
    double argRSGDGradTol = -1;
    double argRSGDRelax = 0.5;

    // parse configure string:
    this->CallArgsInitialize(config, &cla);
    cla.AddArgument("min", this->A_M_SPAC, &argRSGDMin, "");
    cla.AddArgument("scales", this->A_M_SPAC, &argRSGDScales, "");
    cla.AddArgument("steps", this->A_M_SPAC, &argRSGDSteps, "");
    cla.AddArgument("iterations", this->A_M_SPAC, &argRSGDIterations, "");
    cla.AddArgument("gradtol", this->A_SPAC, &argRSGDGradTol, "");
    cla.AddArgument("relax", this->A_SPAC, &argRSGDRelax, "");
    if (!cla.Parse())
      return false; // bad arguments or wrong number of arguments ...

    // cast!
    rsgdOpt = itk::SmartPointer<RSGDOptimizerType>(static_cast<RSGDOptimizerType *>(
      this->m_Optimizer.GetPointer()));

    if (argRSGDMin != 0 && argRSGDMin != 1)
      return false;
    rsgdOpt->SetMinimize(argRSGDMin);

    RSGDOptimizerType::ScalesType weights;
    weights.SetSize(this->m_Transformation->GetNumberOfParameters());
    weights.Fill(1.0); // default
    if (argRSGDScales.size() == this->m_Transformation->GetNumberOfParameters())
    {
      for (unsigned int i = 0; i < argRSGDScales.size(); ++i)
        weights[i] = argRSGDScales[i];
    }
    else if (argRSGDScales.size() != 0)
      return false; // wrong number of arguments
    rsgdOpt->SetScales(weights);

    double minStep = 0.01; // defaults
    double maxStep = 1.0;
    if (argRSGDSteps.size() >= 2)
    {
      minStep = argRSGDSteps[0];
      maxStep = argRSGDSteps[1];
    }
    rsgdOpt->SetMinimumStepLength(minStep);
    rsgdOpt->SetMaximumStepLength(maxStep);

    rsgdOpt->SetNumberOfIterations(argRSGDIterations);

    if (argRSGDGradTol > 0.)
      rsgdOpt->SetGradientMagnitudeTolerance(argRSGDGradTol);

    if (argRSGDRelax <= 0. || argRSGDRelax >= 1.)
      return false;
    rsgdOpt->SetRelaxationFactor(argRSGDRelax);

    return true;
  }
  else
    return Superclass::ConfigureOptimizer(config); // delegate
}

template <typename TInternalPixelType, typename TScalarType>
typename MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType,
  TScalarType>::HistImagePointer
  MultiResolutionImage2D3DRegistrationMethodWithSRC<TInternalPixelType, TScalarType>
::GetCurrentMetricHistogram(bool update)
{
  /* typedefs for concrete casting */
  typedef StochasticRankCorrelationImageToImageMetric<InternalImageType,
       InternalImageType, RankImageType> SRCMetricType;
  typedef typename SRCMetricType::Pointer SRCMetricPointer;

  SRCMetricPointer src = SRCMetricType::New();
  if (this->m_Metric && (this->m_Metric->GetNameOfClass() == src->GetNameOfClass()))
  {
    SRCMetricPointer metric = itk::SmartPointer<SRCMetricType>(
      static_cast<SRCMetricType *>(this->m_Metric.GetPointer()));
    if (update)
      metric->GetValue(metric->GetTransform()->GetParameters());
    // NOTE: metric->GetSampleDistribution() will return the sample distrib.!
    return NULL;
  }
  else
  {
    return Superclass::GetCurrentMetricHistogram(update); // delegate
  }
}


}

#endif ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_TXX_
