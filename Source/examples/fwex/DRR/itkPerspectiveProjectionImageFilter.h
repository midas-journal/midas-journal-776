
#ifndef __ITKPERSPECTIVEPROJECTIONIMAGEFILTER_H
#define __ITKPERSPECTIVEPROJECTIONIMAGEFILTER_H

#include "itkNumericTraits.h"
#include "itkConceptChecking.h"
#include "itkMatrixOffsetTransformBase.h"
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "vnl/vnl_math.h"
#include <vector>
#include <map>
#include <numeric>
#include <utility>

namespace Collector
{

/** \class Collector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Interface for all collectors. Collectors calculate one value \f$v\f$ from many
 * values \f$v_{i}\f$ that get added.
 * These classes are similar to those used in itk::ProjectionImageFilter
 * by Emilian Beronich and Gaetan Lehmann but have been changed to represent
 * "intelligent" pixels of an image.
 *
 * \see Collector::MaximumCollector
 * \see Collector::MinimumCollector
 * \see Collector::SumCollector
 * \see Collector::MeanCollector
 * \see Collector::MedianCollector
 * \see Collector::StandardDeviationCollector
 */
template <class TValue>
class Collector
{
public:
  Collector() {}
  ~Collector() {}

  typedef TValue ValueType;

  /** Initialises and resets all collector variables. */
  virtual inline void Initialize() = 0;

  /** Adds one \a value to the collector.
   * \param value Value is added to the collector.
   */
  virtual inline void Add(const TValue &value) = 0;

  /** Computes the result of the collector based on all added values.
   * \return Computed result of the collector.
   */
  virtual inline TValue GetValue() = 0;
};

/** \class MaximumCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the maximum \f$v_{max}\f$ of all added values \f$v_{i}\f$:
 * \f$v_{max}=\max\left(v_{1},\dots,v_{n}\right)\f$
 */
template <class TValue = double>
class MaximumCollector :
  public Collector<TValue>
{
public:
  MaximumCollector() {}
  ~MaximumCollector() {}

  virtual inline void Initialize()
    {
    m_Maximum = itk::NumericTraits<TValue>::NonpositiveMin();
    }

  virtual inline void Add(const TValue &value)
    {
    m_Maximum = vnl_math_max(m_Maximum, value);
    }

  virtual inline TValue GetValue()
    {
    return m_Maximum;
    }

private:
  TValue m_Maximum;
};

/** \class MinimumCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the minimum \f$v_{min}\f$ of all added values \f$v_{i}\f$:
 * \f$v_{min}=\min\left(v_{1},\dots,v_{n}\right)\f$
 */
template <class TValue = double>
class MinimumCollector :
  public Collector<TValue>
{
public:
  MinimumCollector() {}
  ~MinimumCollector() {}

  virtual inline void Initialize()
    {
    m_Minimum = itk::NumericTraits<TValue>::max();
    }

  virtual inline void Add(const TValue &value)
    {
    m_Minimum = vnl_math_min(m_Minimum, value);
    }

  virtual inline TValue GetValue()
    {
    return m_Minimum;
    }

private:
  TValue m_Minimum;
};

/** \class SumCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the sum \f$v_{sum}\f$ of all added values \f$v_{i}\f$:
 * \f$v_{sum}=\sum_{i=1}^{n}v_{i}\f$
 */
template <class TValue = double>
class SumCollector :
  public Collector<TValue>
{
public:
  SumCollector() {}
  ~SumCollector() {}

  virtual inline void Initialize()
    {
    m_Sum = itk::NumericTraits<TValue>::Zero;
    }

  virtual inline void Add(const TValue &value)
    {
    m_Sum += value;
    }

  virtual inline TValue GetValue()
    {
    return m_Sum;
    }

private:
  TValue m_Sum;
};

/** \class MeanCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the mean \f$\bar{v}\f$ of all added values \f$v_{i}\f$:
 * \f$\bar{v}=\frac{1}{n}\sum_{i=1}^{n}v_{i}\f$
 */
template <class TValue = double, class TCounter = unsigned long>
class MeanCollector :
  public Collector<TValue>
{
public:
  typedef typename itk::NumericTraits<TValue>::RealType TValueRealType;

  MeanCollector() {}
  ~MeanCollector() {}

  virtual inline void Initialize()
    {
    m_Sum = itk::NumericTraits<TValue>::Zero;
    m_Count = itk::NumericTraits<TCounter>::Zero;
    }

  virtual inline void Add(const TValue &value)
    {
    m_Sum += value;
    ++m_Count;
    }

  virtual inline TValueRealType GetValue()
    {
    if (m_Count <= itk::NumericTraits<TCounter>::Zero)
      {
      return itk::NumericTraits<TValue>::Zero;
      }
    else
      {
      return static_cast<TValueRealType>(m_Sum)/m_Count;
      }
    }

private:
  TValue m_Sum;
  TCounter m_Count;
};

/** \class MedianCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the median \f$\tilde{v}\f$ of all added values \f$v_{i}\f$, which
 * are sorted ascending:
 * \f$\tilde{v}=\begin{cases}
 * v_{\frac{n+1}{2}}
 * & \textrm{if }n\textrm{ is odd}\\
 * \frac{1}{2}\left(v_{\frac{n}{2}}+v_{\frac{n}{2}+1}\right)
 * & \textrm{if }n\textrm{ is even}
 * \end{cases}\f$
 */
template <class TValue = double>
class MedianCollector :
  public Collector<TValue>
{
public:
  MedianCollector() {}
  ~MedianCollector() {}

  virtual inline void Initialize()
    {
    m_Values.clear();
    }

  virtual inline void Add(const TValue &value)
    {
    m_Values.push_back(value);
    }

  virtual inline TValue GetValue()
    {
    int size = (int) m_Values.size();
    std::cout <<"s:"<<size<<" ";
    for (int i=0; i<size; i++) std::cout << m_Values[i]<<", ";


    if(size <= 0)
      {
      return  itk::NumericTraits<TValue>::Zero;
      }
    else
      {
      // size is odd
      if (size%2 == 1)
        {
        std::cout <<" odd ";
        typename std::vector<TValue>::iterator medianIt = m_Values.begin() + (size+1)/2 - 1;
        std::nth_element(m_Values.begin(), medianIt, m_Values.end());
        return *medianIt;
        }
      // size is even
      else
        {
        std::cout <<" even ";
        typename std::vector<TValue>::iterator medianIt1 = m_Values.begin() + size/2 -1;
        typename std::vector<TValue>::iterator medianIt2 = medianIt1 + 1;
        std::nth_element(m_Values.begin(), medianIt1, m_Values.end());
        std::nth_element(m_Values.begin(), medianIt2, m_Values.end());
        TValue median = (*medianIt1 + *medianIt2)/2;
        return median;
        }
      }
    }

private:
  std::vector<TValue> m_Values;
};

/** \class StandardDeviationCollector
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * Computes the standard deviation \f$s_{v}\f$ of all added values \f$v_{i}\f$:
 * \f$s_{v}=\sqrt{\frac{1}{n-1}\sum_{i=1}^{n}\left(v_{i}-\bar{v}\right)^{2}}\f
 */
template <class TValue = double>
class StandardDeviationCollector :
  public Collector<TValue>
{
public:
  typedef typename itk::NumericTraits<TValue>::RealType TValueRealType;

  StandardDeviationCollector() {}
  ~StandardDeviationCollector() {}

  virtual inline void Initialize()
    {
    m_Values.clear();
    }

  virtual inline void Add(const TValue &value)
    {
    m_Values.push_back(value);
    }

  virtual inline TValueRealType GetValue()
    {
    typename std::vector<TValue>::size_type count = m_Values.size();
    if(count <= 1)
      {
      return itk::NumericTraits<TValueRealType>::Zero;
      }
    else
      {
      TValue sum = std::accumulate(m_Values.begin(), m_Values.end(),
        itk::NumericTraits<TValue>::Zero);
      TValueRealType mean = static_cast<TValueRealType>(sum)/count;
      TValueRealType variance = itk::NumericTraits<TValueRealType>::Zero;
      typename std::vector<TValue>::iterator it;
      for(it = m_Values.begin(); it != m_Values.end(); it++)
        {
        variance += vnl_math_sqr(*it - mean);
        }
      variance = variance/(count-1);
      return vcl_sqrt(variance);
      }
    }

private:
  std::vector<TValue> m_Values;
};
} // end namespace Collector

namespace itk
{

/** \class PerspectiveProjectionImageFilter
 * \brief Base class for filters that produce an perspectively projected image
 * (digitally reconstructed radiograph, DRR).
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * PerspectiveProjectionImageFilter is the base class for all filters that produce
 * a perspectively projected image (digitally reconstructed radiograph, DRR)
 * of an input image (volume). The volume is transformed with the specified
 * transform and projected onto the output image (DRR), which is positioned in
 * the x/y plane.
 *
 * The class is templated over the types of the input image, output image,
 * the scalar type of the transform and the collector.
 * The collector is  used to collect the projected voxel intensities and to
 * calculate the DRR pixel values. They are used as pixel types of an internal
 * image in the subclasses and offer different types of compositing
 * (e.g. X-ray, MIP).
 *
 * The following parameters must be set to configure the geometry:
 * - Focal point position.
 * - Origin, size and spacing of the DRR (output image). The defaults are
 *   unit spacing, zero origin and zero size.
 * - Origin of the volume.
 * - Parameters of the used transformation (e.g. translation, rotation, rotation
 *   centre).
 *
 * The DRR is positioned in the xy-plane, therefore the origin of the DRR
 * must have a z-coordinate of zero, while the x-coordinate and y-coordinate
 * can be set arbitrary. The DRR is considered as a volumetric dataset with one
 * slice, which means the size in z-direction must be one.
 *
 * The focus (focal point) can be placed anywhere above the DRR with a positive
 * z-coordinate. The volume is positioned in space through its origin and the
 * used transformation.
 * The transformation will be applied to the volume (all voxel positions). All
 * rotations are carried out around a rotation centre that can be placed
 * anywhere in space and must be set in the used transformation.
 * It should be close to the centre of the volume
 * \f$\left(origin+\frac{1}{2}size*spacing\right)\f$  and lie on the ray that
 * travels through the focus and is normal to the xy-plane (normal beam).
 * This avoids rotations that move the volume outside the pyramid that is
 * formed by the focal point and the DRR.
 * The transform is performed in physical space and is applied to every voxel
 * position.
 *
 * Only voxels with intensities above or equal the set threshold will be
 * projected (pre-classification).
 *
 * The intensity transfer function (ITF) is a transfer function that transforms
 * the voxel values (post-classification).
 * It is an stepwise or linear interpolation of values between specified nodes.
 * The nodes represent support points of the ITF and are pairs of input/output
 * intensities. At least two nodes must be set and intensities that are not
 * between two nodes get mapped to zero.
 *
 * This filter overrides several inherited methods from the parent classes
 * to correctly participate in the ITK data processing pipeline.
 * In particular, this filter overrides
 * - ProcessObject::GenerateOutputInformaton()
 * - ProcessObject::GenerateInputRequestedRegion(),
 *   ImageSource::GenerateInputRequestedRegion(),
 *   ImageToImageFilter::GenerateInputRequestedRegion()
 * - Object::GetMTime()
 * - ProcessObject::EnlargeOutputRequestedRegion()
 *
 * \warning The dimension of the input and output image must be three and only
 * scalar voxel/pixel types are supported.
 *
 * \see ImageToImageFilter
 * \see ImageSource
 * \see ProcessObject
 *
 * \ingroup GeometricTransforms IntensityImageFilters
 */
template <class TInputImage, class TOutputImage, class TScalarType=double,
  class TCollector=Collector::SumCollector<double> >
class ITK_EXPORT PerspectiveProjectionImageFilter:
  public ImageToImageFilter<TInputImage, TOutputImage>
{

public:
  /** Standard class typedefs. */
  typedef PerspectiveProjectionImageFilter              Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;

  typedef TOutputImage                             OutputImageType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;

  /** Image typedefs. */
  typedef typename InputImageType::PixelType       InputPixelType;
  typedef typename InputImageType::PointType       InputPointType;
  typedef typename InputImageType::IndexType       InputIndexType;
  typedef typename InputImageType::SizeType        InputSizeType;
  typedef typename InputImageType::RegionType      InputRegionType;
  typedef typename InputImageType::PointType       InputOriginType;
  typedef typename InputImageType::SpacingType     InputSpacingType;

  typedef typename OutputImageType::PixelType      OutputPixelType;
  typedef typename OutputImageType::PointType      OutputPointType;
  typedef typename OutputImageType::IndexType      OutputIndexType;
  typedef typename OutputImageType::SizeType       OutputSizeType;
  typedef typename OutputImageType::RegionType     OutputRegionType;
  typedef typename OutputImageType::PointType      OutputOriginType;
  typedef typename OutputImageType::SpacingType    OutputSpacingType;

  /** Collector typedefs. */
  typedef TCollector                               CollectorType;
  typedef typename TCollector::ValueType           CollectorValueType;

  /** Transform typedefs. */
  typedef TScalarType                              TransformScalarType;
  typedef MatrixOffsetTransformBase<TScalarType,
      itkGetStaticConstMacro(InputImageDimension),
      itkGetStaticConstMacro(InputImageDimension)> TransformType;
  typedef typename TransformType::Pointer          TransformTypePointer;
  typedef typename TransformType::ConstPointer     TransformTypeConstPointer;

  /** Type of Container with nodes of the intensitytransferfunction (ITF) that
   * represent a mapping of pixelvalues from the input image (post-classification).
   * The nodes represent support points of the ITF. Values are interpolated
   * between these nodes dependent on the set method for interpolation.
   * \see SetIntensityTransferFunctionNodes()
   * \see SetIntensityInterpolation()
   * \see IntensityInterpolationEnumType
   */
  typedef std::map<InputPixelType, InputPixelType> IntensityTransferFunctionNodesType;

  /** Enum type that indicates how the interpolation of intensities between
   * the nodes of the intensitytransferfunction (ITF) is performed.
   * \see IntensityTransferFunctionNodesType
   * \see SetIntensityTransferFunctionNodes()
   * \see SetThreshold()
   */
  typedef enum {
    /** The intensities of the inputimage are not transformed. */
    II_NONE,
    /** The intensities of the input image are transformed by a stepwise
     * interpolation between the ITF nodes \f$n_{i}=\left(v_{i},w_{i}\right)\f$.
     * Input values \f$v\f$ between two nodes \f$n_{i}\f$ and \f$n_{i+1}\f$
     * of the ITF are mapped to the output value \f$w_{i}\f$ of the first node
     * \f$n_{i}\f$.
     * If a pixelvalue is outside the specified ITF (not between two nodes)
     * it gets mapped to zero.
     */
    II_STEPWISE,
    /** The intensities of the input image are transformed by a linear
     * interpolation between the ITF nodes \f$n_{i}=\left(v_{i},w_{i}\right)\f$.
     * Input values \f$v\f$ between two nodes \f$n_{i}\f$ and \f$n_{i+1}\f$ of
     * the ITF are mapped to the linear interpolated output value \f$w\f$.
     * The type of the used collector must be \c float or \c double otherwise
     * casting errors will occur.
     * If a pixelvalue is outside the specified ITF (not between two nodes)
     * it gets mapped to zero.
     *
     */
    II_LINEAR} IntensityInterpolationEnumType;

  /** Run-time type information (and related methods). */
  itkTypeMacro(PerspectiveProjectionImageFilter, ImageToImageFilter);

  /** Set the coordinate transformation.
   * Set the coordinate transform to use for transformation of the input image.
   * Note that the transform is performed in the physical space.
   * By default the filter uses an MatrixOffsetTransformBase that is set to
   * identity transform. You must provide a different transform here, before
   * attempting to run the filter.
   */
  itkSetObjectMacro(Transform, TransformType);

  /** Get a pointer to the coordinate transform.
   * \see SetTransform()
   */
  itkGetObjectMacro(Transform, TransformType);

  /** Get a pointer to the coordinate transform.
   * \see SetTransform()
   */
  itkGetConstObjectMacro(Transform, TransformType);

  /** Set the size of the output image (DRR). */
  itkSetMacro(OutputSize, OutputSizeType);

  /** Get the size of the output image (DRR). */
  itkGetConstReferenceMacro(OutputSize, OutputSizeType);

  /** Set the origin of the output image (DRR).
   * Must be in x/y-plane! Therefore the z-coordinate must be zero.
   */
  virtual void SetOutputOrigin(const OutputOriginType);
  /** \copydoc SetOutputOrigin(OutputOriginType) */
  virtual void SetOutputOrigin(const double *);

  /** Get the origin of the output image (DRR). */
  itkGetConstReferenceMacro(OutputOrigin, OutputOriginType);

  /** Set the spacing of the output image (DRR). */
  itkSetMacro(OutputSpacing, OutputSpacingType);
  /** \copydoc SetOutputSpacing(OutputSpacingType) */
  virtual void SetOutputSpacing(const double *);

  /** Get the spacing of the output image (DRR). */
  itkGetConstReferenceMacro(OutputSpacing, OutputSpacingType);

  /** Set focalpoint position for projection. */
  itkSetMacro(FocalPoint, InputPointType);

  /** Get focalpoint position for projection. */
  itkGetConstReferenceMacro(FocalPoint, InputPointType);

  /** Set the threshold. Pixels with intensities above or equal are projected. */
  itkSetMacro(Threshold, InputPixelType);

  /** Get the threshold. */
  itkGetConstMacro(Threshold, InputPixelType);

  /** Set container with nodes of the intensitytransferfunction (ITF).
   * At least two nodes must be in the container otherwise no interpolation
   * is possible.
   * \see SetIntensityTransferFunctionNodesType
   * \see SetIntensityInterpolation()
   * \see IntensityInterpolationEnumType
   */
  virtual void SetIntensityTransferFunctionNodes(const IntensityTransferFunctionNodesType);

  /** Get container with nodes of the intensitytransferfunction (ITF).
   * \see SetIntensityTransferFunctionNodesType
   * \see SetIntensityTransferFunctionNodes()
   */
  virtual const IntensityTransferFunctionNodesType& GetIntensityTransferFunctionNodes(void) const;

  /** Set how the interpolation of intensities between the nodes of the
   * intensitytransferfunction is calculated.
   * \see IntensityInterpolationEnumType
   */
  itkSetEnumMacro(IntensityInterpolation, IntensityInterpolationEnumType);

  /** Get how the interpolation of intensities between the nodes of the
   * intensitytransferfunction is calculated.
   * \see IntensityInterpolationEnumType
   */
  itkGetEnumMacro(IntensityInterpolation, IntensityInterpolationEnumType);

  /** Generate the information describing the output data. This is overridden
   * since the filter produces an image with a different size than its input image.
   * \see ProcessObject::GenerateOutputInformaton()
   */
  virtual void GenerateOutputInformation(void);

  /** Generate the information describing the requested input region in the
   * input image. This is overridden since the filter needs a different input
   * requested region than the output requested region. The entire input image
   * is needed to produce the output.
   * \see ImageToImageFilter::GenerateInputRequestedRegion()
   * \see ImageSource::GenerateInputRequestedRegion()
   * \see ProcessObject::GenerateInputRequestedRegion()
   * */
  virtual void GenerateInputRequestedRegion(void);

  /** Return this objects modified time based on the components. This is
   * overridden since the time when the object was last modified must consider
   * the stored transformation.
   * \see Object::GetMTime()
   */
  unsigned long GetMTime(void) const;

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputPixelType>));
  itkConceptMacro(OutputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<OutputPixelType>));
  /** Input and output images must be the same dimension. */
  itkConceptMacro(ImageDimensionCheck,
    (Concept::SameDimension<itkGetStaticConstMacro(InputImageDimension),
    itkGetStaticConstMacro(OutputImageDimension)>));
  /** End concept checking */
#endif

  /** Stepwise interpolates a mapped value from intensitytransferfunction (ITF)
   * at \a value. If a pixelvalue is outside the specified ITF
   * (not between two nodes) it gets mapped to zero.
   * \param intensityTransferFunction Function with nodes of a function that are
   *              used to map value. Must have at least two nodes but this is
   *              not checked here.
   * \param value Value to map based on intensitytransferfunction.
   * \see IntensityInterpolationEnumType::II_STEPWISE
   * \see IntensityTransferFunctionNodesType
   */
  inline double StepwiseInterpolatePixelValue(
    const IntensityTransferFunctionNodesType &intensityTransferFunction,
    const double &value) const
    {
    double mappedValue = NumericTraits<InputPixelType>::Zero;

    // define iterators to walk intensityTransferFunction (ITF)
    typedef typename IntensityTransferFunctionNodesType::const_iterator ITFIteratorType;
    ITFIteratorType nodeOneIt = intensityTransferFunction.begin();
    ITFIteratorType nodeTwoIt = intensityTransferFunction.begin();
    ++nodeTwoIt;

    // iterate over ITF to find value
    while (nodeTwoIt != intensityTransferFunction.end())
      {
      const double x1 = nodeOneIt->first;
      const double x2 = nodeTwoIt->first;
      if (value >= x1 && value < x2)
        {
        // if value is between two nodes stepwise interpolate and break
        mappedValue = nodeOneIt->second;
        break;
        }
      ++nodeOneIt;
      ++nodeTwoIt;
      }

    return mappedValue;
    }

  /** Linear interpolates a mapped value from intensitytransferfunction (ITF)
   * at \a value. If a pixelvalue is outside the specified ITF
   * (not between two nodes) it gets mapped to zero.
   * \param intensityTransferFunction Function with nodes of a function that are
   *              used to map value. Must have at least two nodes but this is
   *              not checked here.
   * \param value Value to map based on intensitytransferfunction.
   * \see IntensityInterpolationEnumType::II_LINEAR
   * \see IntensityTransferFunctionNodesType
   */
  inline double LinearInterpolatePixelValue(
    const IntensityTransferFunctionNodesType &intensityTransferFunction,
    const double &value) const
    {
    double mappedValue = NumericTraits<InputPixelType>::Zero;

    // define iterators to walk intensityTransferFunction (ITF)
    typedef typename IntensityTransferFunctionNodesType::const_iterator ITFIteratorType;
    ITFIteratorType nodeOneIt = intensityTransferFunction.begin();
    ITFIteratorType nodeTwoIt = intensityTransferFunction.begin();
    ++nodeTwoIt;

    // iterate over ITF to find value
    while (nodeTwoIt != intensityTransferFunction.end())
      {
      const double x1 = nodeOneIt->first;
      const double x2 = nodeTwoIt->first;
      if (value >= x1 && value < x2)
        {
        // if value is between two nodes linear interpolate and break
        const double y1 = nodeOneIt->second;
        const double y2 = nodeTwoIt->second;
        mappedValue = y1 + (y2 - y1) * ((value - x1) / (x2 - x1));
        break;
        }
      ++nodeOneIt;
      ++nodeTwoIt;
      }

    return mappedValue;
    }

protected:
  /** Constructor: Initialise new instance. */
 PerspectiveProjectionImageFilter(void);
  ~PerspectiveProjectionImageFilter(void) {};

  /** Print out a description of self. */
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Indicate that the filter will produce more output than it was requested
   * to produce. This is overridden since the filter needs to produce all of
   * its output. It is not possible to compute only a small region of the DRR,
   * because the whole DRR changes if, for instance, the transform is modified.
   * The output requested region is enlarged to the largest possible region.
   * \see ProcessObject::EnlargeOutputRequestedRegion()
   */
  void EnlargeOutputRequestedRegion(DataObject *data);

private:
  PerspectiveProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Coordinate transform to use.
   * \see SetTransform()
   */
  TransformTypePointer m_Transform;

  /** Size of the output image (DRR). */
  OutputSizeType m_OutputSize;

  /** Origin of the output image (DRR).
   * \see SetOutputOrigin()
   */
  OutputOriginType m_OutputOrigin;

  /** Spacing of the output image (DRR). */
  OutputSpacingType m_OutputSpacing;

  /** Focalpoint position for projection. */
  InputPointType m_FocalPoint;

  /** Pixels with intensities above or equal are projected. */
  InputPixelType m_Threshold;

  /** Get container with nodes of the intensitytransferfunction (ITF).
   * \see SetIntensityTransferFunctionNodesType
   * \see SetIntensityTransferFunctionNodes()
   * \see m_IntensityInterpolation
   */
  IntensityTransferFunctionNodesType m_IntensityTransferFunctionNodes;

  /** Indicates how the interpolation of intensities between the nodes of the
   * intensitytransferfunction are calculated.
   * \see IntensityInterpolationEnumType
   */
  IntensityInterpolationEnumType m_IntensityInterpolation;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkPerspectiveProjectionImageFilter.txx"
#endif

#endif

