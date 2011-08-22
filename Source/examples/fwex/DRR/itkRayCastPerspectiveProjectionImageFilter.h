
#ifndef __ITKRAYCASTPERSPECTIVEPROJECTIONIMAGEFILTER_H
#define __ITKRAYCASTPERSPECTIVEPROJECTIONIMAGEFILTER_H

#include "itkPerspectiveProjectionImageFilter.h"

#include <itkInterpolateImageFunction.h>

namespace itk
{

/** \class RayCastPerspectiveProjectionImageFilter
 * \brief Generates a perspectively projected image (digitally reconstructed
 * radiograph, DRR) with ray casting.
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * RayCastPerspectiveProjectionImageFilter generates a perspectively projected
 * image (DRR) with ray casting.
 * For every pixel of the DRR a ray is cast through the focal point, which
 * travels through the input image (volume). The ray is only cast if it
 * intersects with the volume.
 * In specified distances (step size, sampling distance) values are interpolated
 * with the set interpolation method and added to the current processed pixel in
 * DRR.
 *
 * The inverse transform is used to map the sample locations on the ray into the
 * input image space, therefore the transform is not applied to every voxel
 * position. The transformation is used to map the bounding box of the volume
 * into the transformed space where it is used to check if a sample location is
 * inside the volume. Only sample points that are inside the volume are
 * transformed into the input image space where the interpolation is performed.
 *
 * Interpolated sample values above or equal the set threshold will be added
 * to the corresponding DRR pixel or gets transformed by the ITF (post-classification).
 *
 * The pixel values are limited by the minimum and maximum value of the current
 * pixel type (e.g. char: 0-255).
 *
 * This filter is implemented as a multithreaded filter.
 * It provides the following methods for its implementation:
 * - ImageSource::BeforeThreadedGenerateData()
 * - ImageSource::AfterThreadedGenerateData()
 * - ImageToImageFilter::ThreadedGenerateData()
 *
 * \see PerspectiveProjectionImageFilter
 *
 * \ingroup GeometricTransforms IntensityImageFilters Multithreaded
 */
template <class TInputImage, class TOutputImage, class TScalarType=double,
  class TCollector=Collector::SumCollector<double> >
class ITK_EXPORT RayCastPerspectiveProjectionImageFilter:
  public PerspectiveProjectionImageFilter<TInputImage, TOutputImage, TScalarType, TCollector>
{

public:
  /** Standard class typedefs. */
  typedef RayCastPerspectiveProjectionImageFilter       Self;
  typedef PerspectiveProjectionImageFilter<TInputImage,
  		TOutputImage, TScalarType, TCollector>            Superclass;
  typedef SmartPointer<Self>                            Pointer;
  typedef SmartPointer<const Self>                      ConstPointer;

  /** Some convenient typedefs. */
  itkSuperclassTraitMacro(InputImageType);
  itkSuperclassTraitMacro(InputImagePointer);
  itkSuperclassTraitMacro(InputImageConstPointer);

  itkSuperclassTraitMacro(OutputImageType);
  itkSuperclassTraitMacro(OutputImagePointer);
  itkSuperclassTraitMacro(OutputImageConstPointer);

  /** Image typedefs. */
  itkSuperclassTraitMacro(InputPixelType);
  itkSuperclassTraitMacro(InputPointType);
  itkSuperclassTraitMacro(InputIndexType);
  itkSuperclassTraitMacro(InputSizeType);
  itkSuperclassTraitMacro(InputRegionType);
  itkSuperclassTraitMacro(InputOriginType);
  itkSuperclassTraitMacro(InputSpacingType);

  itkSuperclassTraitMacro(OutputPixelType);
  itkSuperclassTraitMacro(OutputPointType);
  itkSuperclassTraitMacro(OutputIndexType);
  itkSuperclassTraitMacro(OutputSizeType);
  itkSuperclassTraitMacro(OutputRegionType);
  itkSuperclassTraitMacro(OutputOriginType);
  itkSuperclassTraitMacro(OutputSpacingType);

  /** Collector typedefs. */
  itkSuperclassTraitMacro(CollectorType);
  itkSuperclassTraitMacro(CollectorValueType);

  /** Transform typedefs. */
  itkSuperclassTraitMacro(TransformScalarType);
  itkSuperclassTraitMacro(TransformType);
  itkSuperclassTraitMacro(TransformTypePointer);
  itkSuperclassTraitMacro(TransformTypeConstPointer);

  /** Intensitytransferfunction typedefs. */
  itkSuperclassTraitMacro(IntensityTransferFunctionNodesType);
  itkSuperclassTraitMacro(IntensityInterpolationEnumType);

  /** Interpolator typedefs. */
  typedef InterpolateImageFunction<InputImageType, TransformScalarType> InterpolatorType;
  typedef typename InterpolatorType::Pointer       InterpolatorTypePointer;
  typedef typename InterpolatorType::ConstPointer  InterpolatorTypeConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(RayCastPerspectiveProjectionImageFilter, PerspectiveProjectionImageFilter);

  /** Set the interpolator that is used to interpolate values at sample points
   * from the input image (volume).
   * The default is itk::LinearInterpolateImageFunction.
   */
  itkSetObjectMacro(Interpolator, InterpolatorType);

  /** Get a pointer to the interpolator function. */
  itkGetConstObjectMacro(Interpolator, InterpolatorType);

  /** Set step size (sampling distance) between samples on ray. */
  itkSetMacro(RayStepSize, TransformScalarType);

  /** Get step size (sampling distance) between samples on ray. */
  itkGetConstMacro(RayStepSize, TransformScalarType);

  /** This method is used to set the state of the filter before multithreading.
   *
   * All pixel values of the output image (DRR) are set to zero. The size, origin
   * and spacing of the internal DRR are set to the same as the DRR. Each pixel
   * (collector) of the internal DRR is initialised.
   *
   * An Exception is thrown if:
   * - The transform is not set
   * - The inverse transform is not set
   * - The interpolator is not set
   * - The ITF is enabled and has less than two nodes
   *
   * \see ImageSource::BeforeThreadedGenerateData()
   */
  virtual void BeforeThreadedGenerateData(void);

  /** This method is used to set the state of the filter after multithreading.
   *
   * Each pixel value of the DRR is set to the corresponding pixel value
   * (collector value) of the internal DRR. The internal pixel value is limited
   * by the minimum and maximum value of the current pixel type (e.g. char: 0-255).
   *
   * \see ImageSource::AfterThreadedGenerateData()
   */
  virtual void AfterThreadedGenerateData(void);

  /** Return this objects modified time based on the components.  This is
   * overridden since the time when the object was last modified must consider
   * the stored transformation and interpolator.
   * \see Object::GetMTime()
   */
  unsigned long GetMTime(void) const;

protected:
  /** Constructor: Initialise new instance. */
  RayCastPerspectiveProjectionImageFilter(void);
  ~RayCastPerspectiveProjectionImageFilter(void) {};

  /** Print out a description of self. */
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Type of internal output image with collectors as pixel type. It is used to
   * process projected sample values and offers different types of compositing
   * (e.g. X-ray, MIP).
   * \see Collector::Collector
   */
  typedef Image<TCollector,
      itkGetStaticConstMacro(OutputImageDimension)>      InternalOutputImageType;
  typedef typename InternalOutputImageType::Pointer      InternalOutputImageTypePointer;
  typedef typename InternalOutputImageType::ConstPointer InternalOutputImageConstPointer;

  /** Typedef for continuous index of outputimage/internal outputimage. */
  typedef ContinuousIndex<typename OutputPointType::ValueType,
    itkGetStaticConstMacro(OutputImageDimension)>        OutputContinuousIndexType;

  /** This filter is implemented as a multithreaded filter.
   * Therefore, this implementation provides a ThreadedGenerateData() routine which
   * is called for each processing thread. The output image data is allocated
   * automatically by the superclass prior to calling ThreadedGenerateData().
   * ThreadedGenerateData writes only to the portion of the output image
   * specified by the parameter "outputRegionForThread".
   * \see ImageToImageFilter::ThreadedGenerateData()
   * \see ImageToImageFilter::GenerateData()
   * \see ImageSource::GenerateData()
   * \see ImageSource::SplitRequestedRegion()
   */
  void ThreadedGenerateData(const OutputRegionType& outputRegionForThread, int threadId);

private:
  RayCastPerspectiveProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** The inverse coordinate transformation.
   * Coordinate transform to use for transformation of sampe points on the ray
   * from the transformed into the input image space where interpolation is
   * performed.
   * \see PerspectiveProjectionImageFilter::SetTransform()
   */
  TransformTypePointer m_InverseTransform;

  /** Interpolator used for interpolation of sample values from the input image. */
  InterpolatorTypePointer m_Interpolator;

  /** Step size (sampling distance) between samples on the ray. */
  TransformScalarType m_RayStepSize;

  /** Internal representation of the outputimage with collectors as pixels.
   * \see InternalOutputImageType
   * \see Collector::Collector
   */
  InternalOutputImageTypePointer m_InternalOutputImage;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRayCastPerspectiveProjectionImageFilter.txx"
#endif

#endif
