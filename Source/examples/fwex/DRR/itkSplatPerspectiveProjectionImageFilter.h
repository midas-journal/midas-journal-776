
#ifndef __ITKSPLATPERSPECTIVEPROJECTIONIMAGEFILTER_H
#define __ITKSPLATPERSPECTIVEPROJECTIONIMAGEFILTER_H

#include "itkPerspectiveProjectionImageFilter.h"
#include "itkImageRegionIterator.h"

namespace itk
{

/** \class SplatPerspectiveProjectionImageFilter
 * \brief Generates a perspectively projected image (digitally reconstructed
 * radiograph, DRR) with wobbled splatting and focus wobbling.
 * \author Markus Neuner, neuner.markus e_T gmx.net
 * \author phil (modifications), phil.steininger e_T gmail.com
 *
 * SplatPerspectiveProjectionImageFilter takes the voxels from the input image
 * (volume) and transforms their position with the set transform.
 * The transformed voxel positions are projected onto the output image (DRR).
 *
 * The value of the projected voxel is added to the DRR regarding to the set
 * splat-assignement method (nearest, linear).
 *
 * Voxels are only projected if their transformed position is inside
 * the pyramid that is formed by the DRR and the focal point.
 *
 * To suppress artifacts "wobbling" of the focal point with gaussian distributed
 * random offsets is used.
 * How the Gaussian distributed translation of the focus is generated (on the
 * fly or precomputed) depends on the set method (none, pre-continuous,
 * pre-random, random).
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
class ITK_EXPORT SplatPerspectiveProjectionImageFilter:
  public PerspectiveProjectionImageFilter<TInputImage, TOutputImage, TScalarType, TCollector>
{

public:
  /** Standard class typedefs. */
  typedef SplatPerspectiveProjectionImageFilter         Self;
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

  /** Enum type that indicates how the Gaussian distributed translation of
   * the focus is generated (on the fly or precomputed).
   */
  typedef enum {
    /** No focaloffset in x or y direction is performed. The focus position
     * remains constant.
     */
    FOC_NONE,
    /** A specified number of Gaussian distributed random offsets is precomputed
     * and stored in a container. The focal offsets in x and y direction are
     * read in continuous sequential order from this container and starts at the
     * beginning if the end is reached.
     * \see SetFocalRandomOffsetSigma()
     * \see SetFocalRandomOffsetMax()
     * \see SetFocalRandomOffsetQuantity()
     */
    FOC_PRE_CONTINUOUS,
    /** A specified number of Gaussian distributed random offsets is precomputed
     * and stored in a container. The focal offsets in x and y direction are
     * read in random order from this container.
     * \see SetFocalRandomOffsetSigma()
     * \see SetFocalRandomOffsetMax()
     * \see SetFocalRandomOffsetQuantity()
     */
    FOC_PRE_RANDOM,
    /** The Gaussian distributed random offsets in x and y direction are
     * generated "on the fly". No pre-computation is performed.
     * \see SetFocalRandomOffsetSigma()
     * \see SetFocalRandomOffsetMax()
     */
    FOC_RANDOM} FocalOffsetComputationEnumType;

  /** Enum type that specifies how the intensities of projected voxels of the
   * volume are assigned to the pixels of the DRR (splat-assignment). This is
   * needed because generally voxels get projected to off grid locations in the
   * DRR, which requires an inverse interpolation.
   */
  typedef enum {
    /** Assigns the intensity of a projected voxel to the nearest pixel in the
     * DRR (Voronoi Region).
     */
    SA_NEAREST,
    /** Assigns the intensity of a projected voxel to the four surrounding
     * pixels in the DRR (Delaunay Region). The intensity value is distributed
     * based on the distance to the neighbours. The pixel with the smallest
     * distance gets the highest, the one with the greatest distance the
     * smallest part of the voxel intensity. The sum of all four parts is the
     * original intensity value. For instance if the intensity value of one is
     * projected in the middle of four pixels (continuous index of 0.5, 0.5)
     * then 0.25 is added to all four neighbouring pixels.
     *
     * Note that the value type of the used collector must be \c float or
     * \c double,  otherwise errors through casting will occur.
     */
    SA_LINEAR } SplatAssignmentEnumType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(SplatPerspectiveProjectionImageFilter, PerspectiveProjectionImageFilter);

  /** Set the enum type that indicates how the focaloffset computation is performed.
   * \see FocalOffsetComputationEnumType
   */
  itkSetEnumMacro(FocalOffsetComputation, FocalOffsetComputationEnumType);

  /** Get the enum type that indicates how the focaloffset computation is performed.
   * \see FocalOffsetComputationEnumType
   */
  itkGetEnumMacro(FocalOffsetComputation, FocalOffsetComputationEnumType);

  /** Set the standard deviation  (in physical coordinate units) of the Gaussian
   * distributed random translation of the focus parallel to the image plane
   * (xy-plane, DRR).
   * \see SetFocalRandomOffsetMax()
   */
  itkSetMacro(FocalRandomOffsetSigma, double);

  /** Get the standard deviation  (in physical coordinate units) of the Gaussian
   * distributed random translation of the focus parallel to the image plane
   * (xy-plane, DRR).
   * \see SetFocalRandomOffsetSigma()
   */
  itkGetConstMacro(FocalRandomOffsetSigma, double);

  /** Set the number that specifies how many Gaussian distributed random numbers
   * for the translation of the focus are generated if they are pre-computed.
   * Only used at the pre-computation methods pre-continuous and pre-random.
   */
  itkSetMacro(FocalRandomOffsetQuantity, unsigned int);

  /** Get the number that specifies how many Gaussian distributed random numbers
   * for the translation of the focus are generated if they are pre-computed.
   * \see SetFocalRandomOffsetQuantity()
   */
  itkGetConstMacro(FocalRandomOffsetQuantity, unsigned int);

  /** Set if the used random number generator for Gaussian distributed
   * random values is initialised deterministically to produce the same series
   * of random numbers every time a projection is performed. This enables the
   * reproducibility of projections and results.
   * \see FocalOffsetComputationEnumType
   */
  itkSetMacro(ComputeFocalOffsetDeterministically, bool);

  /** Get if the used random number generator for Gaussian distributed
   * random values is initialised deterministically.
   * \see FocalOffsetComputationEnumType
   */
  itkGetConstMacro(ComputeFocalOffsetDeterministically, bool);

  // ComputeFocalOffsetDeterministicallyOn() and ComputeFocalOffsetDeterministicallyOff()
  itkBooleanMacro(ComputeFocalOffsetDeterministically);

  /** Set the enum type that indicates how the intensities of projected pixels
   * are assigned to the pixels of the output image.
   * \see SplatAssignmentEnumType
   */
  itkSetEnumMacro(SplatAssignment, SplatAssignmentEnumType);

  /** Get the enum type that indicates how the intensities of projected pixels
   * are assigned to the pixels of the output image.
   * \see SplatAssignmentEnumType
   */
  itkGetEnumMacro(SplatAssignment, SplatAssignmentEnumType);

  /** This method is used to set the state of the filter before multithreading.
   * \see ImageSource::BeforeThreadedGenerateData()
   */
  virtual void BeforeThreadedGenerateData(void);

  /** This method is used to set the state of the filter after multithreading.
   * \see ImageSource::AfterThreadedGenerateData()
   */
  virtual void AfterThreadedGenerateData(void);

protected:
  /** Constructor: Initialise new instance. */
  SplatPerspectiveProjectionImageFilter(void);
  ~SplatPerspectiveProjectionImageFilter(void) {};

  /** Print out a description of self. */
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** Type of internal output image with collectors as pixel type. It is used to
   * process projected voxel values and offers different types of compositing
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

  /** Type of container for each thread, which saves voxels that get projected
   * outside of the corresponding thread region (\c outputRegionForThread).
   * This can happen because of the Gaussian distributed random translation of
   * the focus. If the projected voxel is inside the DRR but outside the current
   * thread region it must be saved and processed after multi-threading is
   * finished. This must be done because each thread is only allowed to write to
   * its assigned image region.
   *
   * The projected voxel value is added to the according container of the thread
   * with its index in the DRR. When multi-threading is finished
   * (ThreadedGenerateData()) all values are added to the internal DRR
   * (AfterThreadedGenerateData()).
   */
  typedef std::vector<std::pair<OutputIndexType, CollectorValueType> > ThreadOutsideRegionContainerType;

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
  SplatPerspectiveProjectionImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  /** Indicates how the focaloffset computation is performed.
   * \see FocalOffsetComputationEnumType
   */
  FocalOffsetComputationEnumType m_FocalOffsetComputation;

  /** Standard deviation  (in physical coordinate units) of the Gaussian
   * distributed random translation of the focus parallel to the image plane
   * (xy-plane, DRR).
   * \see m_FocalRandomOffsetMax
   */
  double m_FocalRandomOffsetSigma;

  /** Specifies how many Gaussian distributed random numbers for the translation
   * of the focus are generated if they are pre-computed.
   * Only used at the pre-computation methods pre-continuous and pre-random.
   * \see m_FocalOffsetComputation
   * \see m_FocalRandomOffsets
   */
  unsigned int m_FocalRandomOffsetQuantity;

  /** Determines if the used random number generator for Gaussian distributed
   * random values is initialised deterministically to produce the same series
   * of random numbers every time a projection is performed. This enables the
   * reproducibility of projections and results.
   * \see m_FocalOffsetComputation
   */
  bool m_ComputeFocalOffsetDeterministically;

  /** Container that saves the pre-computed Gaussian distributed random values
   * for the translation of the focus.
   * \see m_FocalRandomOffsetQuantity
   * \see m_FocalOffsetComputation
   */
  std::vector<double> m_FocalRandomOffsets;

  /** Indicates how the intensities of projected pixels are assigned to the
   * pixels of the output image.
   * \see SplatAssignmentEnumType
   */
  SplatAssignmentEnumType m_SplatAssignment;

  /** Internal representation of the outputimage with collectors as pixels.
   * \see InternalOutputImageType
   * \see Collector::Collector
   */
  InternalOutputImageTypePointer m_InternalOutputImage;

  /** Type of container that is used to store pixels that get projected outside
   * of the outputRegionForThread because of the focal offset in x or y direction.
   * Container that saves one container for each thread, which saves voxels that
   * get projected outside of the corresponding thread region
   * (\c outputRegionForThread). This can happen because of the Gaussian
   * distributed random translation of the focus.
   * \see ThreadOutsideRegionContainerType
   */
  std::vector<ThreadOutsideRegionContainerType > m_ThreadOutsideRegionPixels;

  /** Adds \a value to \a index through iterator \a it on the internal outputimage.
   * \param index Index in m_InternalOutputImage where value is set.
   * \param value Value to set at index.
   * \param it    Iterator for m_InternalOutputImage to access the specified
   *              index and add value efficiently.
   * \see m_InternalOutputImage
   * \see InternalOutputImageType
   * \see Collector::Collector
   */
  inline void AddToInternalOutputImage(const OutputIndexType &index,
    const double &value, ImageRegionIterator<InternalOutputImageType> &it) const
    {
    it.SetIndex(index);
    it.Value().Add(value);
    return;
    }

  /** Checks if index is inside the outputRegionForThread.
   * If \a index is inside \a value is added to \a index through iterator \a it
   * on the internal outputimage (AddToInternalOutputImage()).
   * If \a index is outside \a index and \a value are added to
   * m_ThreadOutsideRegionPixels to be handled when ThreadedGenerateData()
   * is finished.
   * \param index Index in m_InternalOutputImage where value is set.
   * \param value Value to set at index.
   * \param it    Iterator for m_InternalOutputImage to access the specified
   *              index and add value efficiently.
   * \param outputRegionForThread Region of thread with threadId.
   * \param threadId Number of thread.
   * \see AddToInternalOutputImage()
   * \see ThreadedGenerateData()
   */
  inline void AddToThreadInternalOutputImage(const OutputIndexType &index,
    const double &value, ImageRegionIterator<InternalOutputImageType> &it,
    const OutputRegionType &outputRegionForThread, const int &threadId)
    {
    if (outputRegionForThread.IsInside(index))
      {
      this->AddToInternalOutputImage(index, value , it);
      }
    else
      {
      this->m_ThreadOutsideRegionPixels[threadId].push_back(std::make_pair(index, value));
      }
    return;
    }

  /** Assigns the intensity of a projected voxel to the four surrounding pixels
   * in the DRR (Delaunay Region). The intensity value is distributed based on
   * the distance to the neighbours.
   * \param index Continuous index in \c m_InternalOutputImage where value is projected.
   * \param value Value to set at index.
   * \param it    Iterator for \c m_InternalOutputImage to access the specified
   *              index and add value efficiently.
   * \param outputRegionForThread Region of thread with threadId.
   * \param threadId Number of thread.
   * \see AddToThreadInternalOutputImage()
   * \see SplatAssignmentEnumType::SA_LINEAR
   * \see m_SplatAssignment
   */
  inline void SplatAssignmentLinear(const OutputContinuousIndexType &index,
    const double &value, ImageRegionIterator<InternalOutputImageType> &it,
    const OutputRegionType &outputRegionForThread, const int &threadId)
    {
    typedef typename OutputIndexType::IndexValueType IndexValueType;
    typedef typename OutputContinuousIndexType::ValueType ContinuousIndexValueType;
    const IndexValueType indexValueZero = NumericTraits<IndexValueType>::Zero;

    // The four neighbours are numbered binary with 00=base, 01, 10, 11
    // compute index of neighbour 00
    OutputIndexType indexBase;

    // only x and y dimension is needed
    ContinuousIndexValueType indexDimension = index[0];
    indexBase[0] = static_cast<IndexValueType>(indexDimension);
    if(indexDimension < indexValueZero
      && ContinuousIndexValueType(indexBase[0]) != indexDimension)
      {
      indexBase[0]--;
      }

    indexDimension = index[1];
    indexBase[1] = static_cast<IndexValueType>(indexDimension);
    if(indexDimension < indexValueZero
      && ContinuousIndexValueType(indexBase[1]) != indexDimension)
      {
      indexBase[1]--;
      }

    indexBase[2] = indexValueZero;

    // get distance of projected pixel from base-index in x and y dimension
    const double distanceX = index[0] - double(indexBase[0]);
    const double distanceY = index[1] - double(indexBase[1]);

    // distribute value based on distance to neighbours
    const double inputPixelValueA = value*(1-distanceY);
    const double inputPixelValueB = value*distanceY;

    // set value of neighbour 00
    this->AddToThreadInternalOutputImage(indexBase, inputPixelValueA*(1-distanceX),
        it, outputRegionForThread, threadId);
    // set value of neighbour 10
    ++indexBase[0];
    this->AddToThreadInternalOutputImage(indexBase, inputPixelValueA*distanceX,
        it, outputRegionForThread, threadId);
    // set value of neighbour 11
    ++indexBase[1];
    this->AddToThreadInternalOutputImage(indexBase, inputPixelValueB*distanceX,
        it, outputRegionForThread, threadId);
    // set value of neighbour 01
    --indexBase[0];
    this->AddToThreadInternalOutputImage(indexBase, inputPixelValueB*(1-distanceX),
        it, outputRegionForThread, threadId);
    return;
    }

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSplatPerspectiveProjectionImageFilter.txx"
#endif

#endif
