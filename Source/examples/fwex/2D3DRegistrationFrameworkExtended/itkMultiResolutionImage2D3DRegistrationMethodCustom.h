

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_H_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_H_


#include <itkMultiResolutionImage2D3DRegistrationMethod.h>
#include <itkMultiResolutionRegistrationCommandCustom.h>


namespace itk
{


/** \class MultiResolutionImage2D3DRegistrationMethodCustom
 *  \brief This class extends the 2D/3D-registration framework.
 *
 * <p>This class extends the 2D/3D-registration framework. A new metric type
 * is supported: Mattes Mutual Information (2D/3D-registration compatible
 * version). </p>
 *
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.4
 */
template <typename TInternalPixelType = float, typename TScalarType = double>
class MultiResolutionImage2D3DRegistrationMethodCustom :
  public MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType, 
  TScalarType>
{
public:
  /** Standard class typedefs. */
  typedef MultiResolutionImage2D3DRegistrationMethodCustom Self;
  typedef MultiResolutionImage2D3DRegistrationMethod<TInternalPixelType,
    TScalarType> Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** typedefs from superclass **/
  typedef typename Superclass::HistImageType HistImageType; 
  typedef typename Superclass::HistImagePointer HistImagePointer;
  typedef typename Superclass::InternalImageType InternalImageType;

  /** new multi-resolution command **/
  typedef MultiResolutionRegistrationCommand<Self>
    MultiResolutionObserverType;
  typedef typename MultiResolutionObserverType::Pointer
    MultiResolutionObserverPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(MultiResolutionImage2D3DRegistrationMethodCustom,
    MultiResolutionImage2D3DRegistrationMethod);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  /**
   * Configure the internal metric by applying a configuration string.
   * (This is maybe more comfortable than accessing the components externally.)
   * @param config string expression of format:
   *
   * @see itk::MultiResolutionImage2D3DRegistrationMethod::ConfigureMetric()
   *
   * additionally:
   *
   * (MMI)
   * "histbins <int> lowerbounds <double> <double>
   * upperbounds <double> <double> nosamples <double>"
   *
   * where
   *
   * <b>histbins</b> (integer) specifies the size of the histogram
   * (moving & fixed have the same number of bins),
   * <b>lowerbounds</b> (doubles) specifies the lower intensity bounds of
   * the histogram (moving, fixed),
   * <b>upperbounds</b> (doubles) specifies the upper intensity bounds of the
   * histogram (moving, fixed),
   * <b>nosamples</b> (double) specifies the number of samples (rounded, if > 1) 
   * or the percentage of the sample pixels in fixed image region (0.0..1.0) or 
   * all pixels (if < 0).
   *
   * @result true if configure string has successfully been applied
   **/
  virtual bool ConfigureMetric(const std::string config);

  /**
   * @override itk::MultiResolutionImage2D3DRegistrationMethod::
   * GetCurrentMetricHistogram(bool) in order to support MMI-joint-PDF too
   */
  virtual HistImagePointer GetCurrentMetricHistogram(bool update);

protected:
  /** Hidden constructor. **/
  MultiResolutionImage2D3DRegistrationMethodCustom();

  /** Discarded destructor. **/
  virtual ~MultiResolutionImage2D3DRegistrationMethodCustom();

  /** Print class information. **/
  virtual void PrintSelf(std::ostream& os, Indent indent) const;

private:
  // purposely not implemented:
  MultiResolutionImage2D3DRegistrationMethodCustom(const Self&); 
  // purposely not implemented:
  void operator=(const Self&);

};


}


#include "itkMultiResolutionImage2D3DRegistrationMethodCustom.txx"


#endif /* ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCUSTOM_H_ */
