//
#ifndef ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_H_
#define ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_H_

#include <itkMultiResolutionImage2D3DRegistrationMethodCustom.h>

#include "oraMultiResolutionImage2D3DRegistrationMethodCommandWithSRC.h"
#include "oraMultiResolutionRegistrationCommandWithSRC.h"
#include "oraOptimizerIterationCommandWithSRC.h"

namespace ora /** open radART **/
{

/** \class MultiResolutionImage2D3DRegistrationMethodWithSRC
 *  \brief This class extends the 2D/3D registration framework for X-rays.
 *
 * This class extends the 2D/3D registration framework for X-rays by providing
 * additional support for<br>
 * - stochastic rank correlation (SRC) metric, <br>
 * - regular step gradient descent optimization, <br>
 * - and fixed image masks.<br>
 *
 * NOTE:<br>
 * If SRC is configured, the PDF-related keywords are interpreted differently.
 * Instead of a PDF, the "sample distribution" (fixed image samples vs.
 * moving image samples are written to comma-separated data files).
 *
 * @see ora::StochasticRankCorrelationImageToImageMetric
 * @see itk::RegularStepGradientDescentBaseOptimizer
 * @see itk::MultiResolutionImage2D3DRegistrationMethodCustom
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.1
 */
template<typename TInternalPixelType = float, typename TScalarType = double>
class MultiResolutionImage2D3DRegistrationMethodWithSRC :
    public itk::MultiResolutionImage2D3DRegistrationMethodCustom<
        TInternalPixelType, TScalarType>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionImage2D3DRegistrationMethodWithSRC Self;
  typedef itk::MultiResolutionImage2D3DRegistrationMethodCustom<TInternalPixelType,
    TScalarType> Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  typedef itk::SmartPointer<const Self> ConstPointer;

  /** typedefs from superclass **/
  typedef typename Superclass::InternalImageType InternalImageType;
  typedef typename Superclass::HistImagePointer HistImagePointer;
  typedef unsigned int RankPixelType;
  typedef itk::Image<RankPixelType, 3> RankImageType;

  /** observer types **/
  typedef ora::MultiResolutionImage2D3DRegistrationMethodCommandWithSRC<Self>
    GenericFrameworkObserverType;
  typedef typename GenericFrameworkObserverType::Pointer
    GenericFrameworkObserverPointer;
  typedef MultiResolutionRegistrationCommandWithSRC<Self>
    MultiResolutionObserverType;
  typedef typename MultiResolutionObserverType::Pointer
    MultiResolutionObserverPointer;
  typedef OptimizerIterationCommandWithSRC<Self> OptimizationObserverType;
  typedef typename OptimizationObserverType::Pointer
    OptimizationObserverPointer;

   /** New mask image type **/
  typedef itk::Image<unsigned char, InternalImageType::ImageDimension>
    FixedMaskImageType;
  typedef typename FixedMaskImageType::Pointer FixedMaskImagePointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionImage2D3DRegistrationMethodWithSRC,
    MultiResolutionImage2D3DRegistrationMethodCustom);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /**
   * Configure the internal metric by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * @see itk::MultiResolutionImage2D3DRegistrationMethodCustom::ConfigureMetric()
   *
   * additionally:
   *
   * (SRC)
   * "rseeds <int> <int> dscales <double> ... <double>
   * fhist <double> <double> <int> <bool> mhist <double> <double> <int> <bool>
   * coverage <double> horn <bool> zerorankscontribute <bool>"
   *
   * where
   *
   * <b>rseeds</b> (integer) specifies random seeds for deterministic behavior
   * of internal stochastic mask generation (if >0) - by default all rseeds are
   * set to 0 (non-deterministic behavior),
   * <b>dscales</b> (double) is a vector of scales for derivative computation
   * (finite distances per transformation parameter - we need
   * NumberOfTransformationParameters scales!),
   * <b>fhist</b> defines the fixed histogram settings: minimum intensity,
   * maximum intensity, number of bins, clip-at-ends-flag (REQUIRED),
   * <b>mhist</b> defines the moving histogram settings: minimum intensity,
   * maximum intensity, number of bins, clip-at-ends-flag (REQUIRED),
   * <b>coverage</b> (double) specifies the sample coverage (between 0 and 100 -
   * specified in percentage of fixed image region),
   * <b>horn</b> (bool) specifies whether Horn-correction for tied ranks should
   * be used,
   * <b>zerorankscontribute</b> (bool) specifies whether 0-ranks of the moving
   * image contribute to SRC coefficient computation.
   *
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *
   * FURTHERMORE, this method introduces a new additional argument applicable
   * to ALL metrics! The effect of this (fixed) mask will however depend on the
   * metric's capabilities. It is not guaranteed that each metric supports the
   * fixed image mask!
   *
   * "circularfixedmask <int> <int> <int>"
   *
   * where
   *
   * <b>circularfixedmask</b> (floats) defines a circular fixed image mask with
   * the first two ints defining the x and y coordinates of the center, and
   * the last int being the radius (all specifications in mm, center is
   * expressed RELATIVE to the FIXED IMAGE origin).
   *
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureMetric(const std::string config);

  /**
   * Configure the internal optimizer by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * @see itk::MultiResolutionImage2D3DRegistrationMethod::ConfigureOptimizer()
   *
   * (REGSTEPGRADDESC)
   * "min <bool> scales <double> ... <double> steps <double> <double>
   * iterations <int> gradtol <double> relax <double>"
   *
   * where
   *
   * <b>min</b> (bool) specifies the optimization criterion (0=maximize,
   * 1=minimize),
   * <b>scales</b> (doubles) specifies the weightings of the
   * transformation dimensions,
   * <b>steps</b> (doubles) specifies the minimum and maximum step lengths of
   * gradient descent,
   * <b>iterations</b> (integer) specifies the maximum number of iterations to
   * be performed (unless gradient tolerance is undershot),
   * <b>gradtol</b> (double) specifies the gradient magnitude tolerance,
   * <b>relax</b> (double) specifies the step length relaxation factor (between
   * 0.0 and 1.0) which is applied if gradient direction changes.
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureOptimizer(const std::string config);

  /**
   * "Misuse" this mechanism for sample distribution extraction (SRC).
   * @param update if TRUE then metric->GetValue() is called to ensure that the
   * sample distribution corresponds to current transformation settings;
   * otherwise the metric is assumed to be up-to-date!
   * @return always NULL (as we do not extract a histogram here!)
   */
  virtual HistImagePointer GetCurrentMetricHistogram(bool update);

protected:
  /** Hidden constructor. **/
  MultiResolutionImage2D3DRegistrationMethodWithSRC();

  /** Discarded destructor. **/
  virtual ~MultiResolutionImage2D3DRegistrationMethodWithSRC();

  /** Print class information. **/
  virtual void PrintSelf(std::ostream& os, itk::Indent indent) const;

private:
  // purposely not implemented:
  MultiResolutionImage2D3DRegistrationMethodWithSRC(const Self&);
  // purposely not implemented:
  void operator=(const Self&);

};

}

#include "oraMultiResolutionImage2D3DRegistrationMethodWithSRC.txx"

#endif /* ORAMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODWITHSRC_H_ */
