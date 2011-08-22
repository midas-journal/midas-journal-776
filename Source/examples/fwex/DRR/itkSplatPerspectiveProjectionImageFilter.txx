
#ifndef __ITKSPLATPERSPECTIVEPROJECTIONIMAGEFILTER_TXX
#define __ITKSPLATPERSPECTIVEPROJECTIONIMAGEFILTER_TXX

#include "itkSplatPerspectiveProjectionImageFilter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkContinuousIndex.h"
#include "itkRealTimeClock.h"
#include "vnl/vnl_random.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
SplatPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::SplatPerspectiveProjectionImageFilter()
{
  m_FocalOffsetComputation = FOC_PRE_CONTINUOUS;
  m_FocalRandomOffsetSigma = 1.0;
  m_FocalRandomOffsetQuantity = 50;
  m_ComputeFocalOffsetDeterministically = false;
  // m_FocalRandomOffsets see BeforeThreadedGenerateData()
  m_SplatAssignment = SA_NEAREST;
  m_InternalOutputImage = InternalOutputImageType::New();
  // m_ThreadOutsideRegionPixels see BeforeThreadedGenerateData()
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
SplatPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "FocalOffsetComputation: " << m_FocalOffsetComputation << std::endl;
  os << indent << "FocalRandomOffsetSigma: "
     << static_cast<typename NumericTraits<double>::PrintType>(m_FocalRandomOffsetSigma)
     << std::endl;
  os << indent<< "FocalRandomOffsetQuantity: "
     << static_cast<typename NumericTraits<unsigned int>::PrintType>(m_FocalRandomOffsetQuantity)
     << std::endl;
  os << indent<< "ComputeFocalOffsetDeterministically: "
     << static_cast<typename NumericTraits<unsigned int>::PrintType>(m_ComputeFocalOffsetDeterministically)
     << std::endl;
  os << indent << "SplatAssignment: " << m_SplatAssignment << std::endl;
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
SplatPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::BeforeThreadedGenerateData(void)
{
  // Get members from superclass
  TransformTypeConstPointer transform = this->GetTransform();
  const IntensityTransferFunctionNodesType& intensityTransferFunctionNodes =
      this->GetIntensityTransferFunctionNodes();
  const IntensityInterpolationEnumType& intensityInterpolation =
      this->GetIntensityInterpolation();

  if(!transform)
    {
    itkExceptionMacro(<< "Transform not set");
    }

  if(intensityInterpolation != Superclass::II_NONE && intensityTransferFunctionNodes.size() < 2)
    {
    itkExceptionMacro(<< "IntensityTransferFunctionNodes has not enough values to interpolate (<2).");
    }

  // Get pointer to output image and initialise image pixels
  OutputImagePointer outputPtr = this->GetOutput();
  outputPtr->FillBuffer(NumericTraits<OutputPixelType>::Zero);

  // Initialise and set up internal outputimage
  // copy region, origin, spacing from outputimage
  m_InternalOutputImage->SetRegions(outputPtr->GetLargestPossibleRegion());
  m_InternalOutputImage->SetOrigin(outputPtr->GetOrigin());
  m_InternalOutputImage->SetSpacing(outputPtr->GetSpacing());
  m_InternalOutputImage->Allocate();

  // Create iterator to access internal outputimage
  typedef ImageRegionIterator<InternalOutputImageType> InternalIteratorType;
  InternalIteratorType internalIt = InternalIteratorType(m_InternalOutputImage,
      m_InternalOutputImage->GetLargestPossibleRegion());

  // Walk internal outputimage and initialise collectors
  internalIt.GoToBegin();
  while (!internalIt.IsAtEnd())
    {
    internalIt.Value().Initialize();
    ++internalIt;
    }

  // initialise container for projected pixels that lie outside the thread's region
  unsigned int numberOfThreads = this->GetNumberOfThreads();
  if(m_ThreadOutsideRegionPixels.size() != numberOfThreads)
    {
    // initialise container
    m_ThreadOutsideRegionPixels.clear();
    m_ThreadOutsideRegionPixels.reserve(numberOfThreads);

    // each thread gets a container for index/value pairs
    m_ThreadOutsideRegionPixels = std::vector<ThreadOutsideRegionContainerType >
      (numberOfThreads, ThreadOutsideRegionContainerType());

    // Reserve an initial capacity for each thread
    for (typename std::vector<ThreadOutsideRegionContainerType >::iterator it
        = m_ThreadOutsideRegionPixels.begin(); it != m_ThreadOutsideRegionPixels.end(); ++it)
      {
      it->reserve(128);
      }
    }

  // precompute gaussian distributed random offsets for focalpoint
  if((m_FocalOffsetComputation == FOC_PRE_RANDOM
    || m_FocalOffsetComputation == FOC_PRE_CONTINUOUS)
    && m_FocalRandomOffsets.size() != m_FocalRandomOffsetQuantity)
    {
    // random number generator
    vnl_random rand;
    if(m_ComputeFocalOffsetDeterministically)
      {
      rand.reseed(1234567);
      }

    // initialise container
    m_FocalRandomOffsets.clear();
    m_FocalRandomOffsets.reserve(m_FocalRandomOffsetQuantity);

    // fill container with gaussian distributed random-values
    while(m_FocalRandomOffsets.size() < m_FocalRandomOffsetQuantity)
      {
      // scale gaussian distributed randomnumber from N(0,1) to N(0,Sigma)
      const double randomVariable = rand.normal() * m_FocalRandomOffsetSigma;

      // add value to container
      m_FocalRandomOffsets.push_back(randomVariable);
      }
    }

  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
SplatPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::AfterThreadedGenerateData(void)
{
  // Get the output pointer
  OutputImagePointer outputPtr = this->GetOutput();
  unsigned int numberOfThreads = this->GetNumberOfThreads();

  // Create iterator to access internal outputimage
  typedef ImageRegionIterator<InternalOutputImageType> InternalIteratorType;
  typename InternalOutputImageType::RegionType maxReg = m_InternalOutputImage->
    GetLargestPossibleRegion();
  InternalIteratorType internalIt = InternalIteratorType(m_InternalOutputImage,
      maxReg);
  typename InternalOutputImageType::IndexType minInd = maxReg.GetIndex();
  typename InternalOutputImageType::IndexType maxInd;
  for (unsigned int i = 0; i < 3; i++)
    {
    maxInd[i] = minInd[i] + maxReg.GetSize()[i] - 1;
    if (maxInd[i] < minInd[i])
      maxInd[i] = minInd[i];
    }

  // Add all pixels that were projected outside the thread's region
  for(unsigned int i = 0; i < numberOfThreads; i++)
    {
    // Get pointer to container of thread i
    ThreadOutsideRegionContainerType *threadContainer = &m_ThreadOutsideRegionPixels[i];

    // Add all pixels in container to internal outputimage
    while(!threadContainer->empty())
      {
      const OutputIndexType index = threadContainer->back().first;

      // NOTE:
      // there are indices outside range of internal output image's region (do
      // not know why at the moment); so just check the bounds and forget
      // invalid pixels:
      if (index[0] >= minInd[0] && index[0] <= maxInd[0] &&
          index[1] >= minInd[1] && index[1] <= maxInd[1] &&
          index[2] >= minInd[2] && index[2] <= maxInd[2])
        {
        const double value = threadContainer->back().second;
        this->AddToInternalOutputImage(index, value , internalIt);
        }

      threadContainer->pop_back();
      }
    }

  // minimum/maximum values of the OutputPixelType
  const OutputPixelType outputMinPixelValue = NumericTraits<OutputPixelType>::NonpositiveMin();
  const OutputPixelType outputMaxPixelValue = NumericTraits<OutputPixelType>::max();
  const CollectorValueType minPixelValue = static_cast<CollectorValueType>(outputMinPixelValue);
  const CollectorValueType maxPixelValue = static_cast<CollectorValueType>(outputMaxPixelValue);

  // Create iterator to access outputimage
  typedef ImageRegionIterator<OutputImageType> OutputIteratorType;
  OutputIteratorType outputIt = OutputIteratorType(outputPtr,
      outputPtr->GetLargestPossibleRegion());

  // get values from internal outputimage and set them in outputimage
  internalIt.GoToBegin();
  outputIt.GoToBegin();

  while (!internalIt.IsAtEnd())
    {
    // get value from internal outputimage
    const CollectorValueType pixelValue = internalIt.Value().GetValue();

    // check if value is inside OutputPixelType and set it
    if(pixelValue > maxPixelValue)
      {
      outputIt.Set(outputMaxPixelValue);
      }
    else if(pixelValue < minPixelValue)
      {
      outputIt.Set(outputMinPixelValue);
      }
    else
      {
      outputIt.Set(static_cast<OutputPixelType>(pixelValue));
      }
    ++internalIt;
    ++outputIt;
    }
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
SplatPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::ThreadedGenerateData(const OutputRegionType &outputRegionForThread, int threadId)
{
  // Get the input and output ConstPointer
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImageConstPointer outputPtr = this->GetOutput();

  // Get members from superclass
  TransformTypeConstPointer transform = this->GetTransform();
  // const OutputSizeType& outputSize = this->GetOutputSize();
  // const OutputOriginType& outputOrigin = this->GetOutputOrigin();
  const OutputSpacingType& outputSpacing = this->GetOutputSpacing();
  const InputPointType& focalPoint = this->GetFocalPoint();
  const InputPixelType& threshold = this->GetThreshold();
  const IntensityTransferFunctionNodesType& intensityTransferFunctionNodes =
      this->GetIntensityTransferFunctionNodes();
  const IntensityInterpolationEnumType& intensityInterpolation =
      this->GetIntensityInterpolation();

  // Create iterator to access inputimage
  typedef ImageRegionConstIteratorWithIndex<InputImageType> InputIteratorType;
  InputIteratorType inputIt(inputPtr, inputPtr->GetRequestedRegion());

  // Create iterator to access internal outputimage
  typedef ImageRegionIterator<InternalOutputImageType> InternalIteratorType;
  InternalIteratorType internalIt = InternalIteratorType(m_InternalOutputImage,
      outputRegionForThread);

  // Create iterator to access randomnumbers
  std::vector<double>::const_iterator randomIt = m_FocalRandomOffsets.begin();

  // Define a few points that will be used to project an input pixel onto the
  // outputimage that is positioned on the x/y plane

  // Cartesian coordinates of current input pixel
  InputPointType inputPoint;
  inputPoint.Fill(0);

  // Transformed Cartesian coordinates of current input pixel
  typename TransformType::OutputPointType transformedInputPoint;
  transformedInputPoint.Fill(0);

  // Cartesian coordinates of projected output pixel
  OutputPointType outputPoint;
  outputPoint.Fill(0);

  // Continuous index of projected output pixel
  typedef ContinuousIndex<typename OutputPointType::ValueType,
      itkGetStaticConstMacro(OutputImageDimension)> OutputContinuousIndexType;
  OutputContinuousIndexType outputContinuousIndex;
  outputContinuousIndex.Fill(0);

  // Index of projected output pixel
  OutputIndexType outputIndex;
  outputIndex.Fill(0);

  // coordinate components of focalpoint
  const TransformScalarType focalX = static_cast<TransformScalarType>(focalPoint[0]);
  const TransformScalarType focalY = static_cast<TransformScalarType>(focalPoint[1]);
  const TransformScalarType focalZ = static_cast<TransformScalarType>(focalPoint[2]);

  // Define a few points that will be used to check if a transformed point is
  // inside the pyramid that is formed by the region and the focal point

  // Get indices of corners from current thread's region
  typedef typename OutputRegionType::IndexType OutputRegionIndexType;
  typedef typename OutputRegionType::SizeType OutputRegionSizeType;
  typedef typename OutputRegionType::IndexValueType OutputRegionIndexValueType;

  // smallest index
  const OutputRegionIndexType regionIndexMin = outputRegionForThread.GetIndex();
  const OutputRegionSizeType regionSize = outputRegionForThread.GetSize();

  // highest index (only x and y dimension is relevant)
  OutputRegionIndexType regionIndexMax;
  regionIndexMax[0] = regionIndexMin[0] +
      static_cast<OutputRegionIndexValueType>(regionSize[0]) - 1;
  regionIndexMax[1] = regionIndexMin[1] +
      static_cast<OutputRegionIndexValueType>(regionSize[1]) - 1;
  regionIndexMax[2] = 0;

  // Vertices of region for current thread (only x and y dimension is relevant)
  // These points are used to check if point is inside pyramid
  // pyramid faces touch because of +-(m_OutputSpacing/2.0)
  OutputPointType regionPointMin;
  outputPtr->TransformIndexToPhysicalPoint(regionIndexMin, regionPointMin);
  regionPointMin = regionPointMin - (outputSpacing/2.0);
  regionPointMin[2] = 0;

  OutputPointType regionPointMax;
  outputPtr->TransformIndexToPhysicalPoint(regionIndexMax, regionPointMax);
  regionPointMax = regionPointMax + (outputSpacing/2.0);
  regionPointMax[2] = 0;

  // Triangle vertices to test if transformed pixel is inside pyramid
  // pyramid check is based on two triangle inside checks (x/z-plane and y/z-plane)
  // based on linear combination of the corner point vectors

  // Determine if x/z-plane or y/z-plane is checked first
  // The dimension with the greater size is checked first
  bool checkYZPlaneFirst = true;  // y/z-plane first
  if(regionSize[0]<regionSize[1])
  {
    checkYZPlaneFirst = false;  // x/z-plane first
  }

  // Triangle 1 vertex-components
  const double A1x = regionPointMin[checkYZPlaneFirst];
  const double B1x = regionPointMax[checkYZPlaneFirst];
  const double C1x = focalPoint[checkYZPlaneFirst];

  // Triangle 2 vertex-components
  const double A2x = regionPointMin[!checkYZPlaneFirst];
  const double B2x = regionPointMax[!checkYZPlaneFirst];
  const double C2x = focalPoint[!checkYZPlaneFirst];

  // in both triangles the same
  const double C12y = focalPoint[2];

  // Precomputed terms for pyramid-check
  const double invDenominator1 = 1 / (C12y*B1x - A1x*C12y);
  const double invDenominator2 = 1 / (C12y*B2x - A2x*C12y);

  // random number generator (needed when m_FocalOffsetComputation = FOC_RANDOM)
  vnl_random rand;
  if(m_ComputeFocalOffsetDeterministically)
    {
    unsigned long seed = 1234567 + (threadId+1) * (outputRegionForThread.GetNumberOfPixels()+1);
    rand.reseed(seed);
    }

  // Iterate over inputimage and project pixels
  inputIt.GoToBegin();

  while (!inputIt.IsAtEnd())
    {
    // Check if pixelvalue is above threshold
    const InputPixelType inputPixelValue = inputIt.Get();
    if (inputPixelValue >= threshold)
      {
      // Get the coordinates of input pixel and perform the transform
      inputPtr->TransformIndexToPhysicalPoint(inputIt.GetIndex(), inputPoint);
      transformedInputPoint = transform->TransformPoint(inputPoint);

      // Test if point is inside pyramid (triangle1)
      double Px = transformedInputPoint[checkYZPlaneFirst];
      const double Py = transformedInputPoint[2];
      double beta =  (A1x*Py - A1x*C12y + C12y*Px - C1x*Py) * invDenominator1;
      double gamma = (B1x*Py - A1x*Py) * invDenominator1;

      if ((beta > 0) && (gamma > 0) && (beta + gamma < 1))
        {
        // Test if point is inside pyramid (triangle2)
        Px = transformedInputPoint[!checkYZPlaneFirst];
        beta =  (A2x*Py - A2x*C12y + C12y*Px - C2x*Py ) * invDenominator2;
        gamma = (B2x*Py - A2x*Py) * invDenominator2;

        if ((beta > 0) && (gamma > 0) && (beta + gamma < 1))
          {
          // Project transformed point onto x/y-plane
          const TScalarType pointZ = transformedInputPoint[2];
          switch (m_FocalOffsetComputation)
            {
            case FOC_NONE:
            default:
              outputPoint[0] = (pointZ*focalX - transformedInputPoint[0]*focalZ)
                                / (pointZ-focalZ);
              outputPoint[1] = (pointZ*focalY - transformedInputPoint[1]*focalZ)
                                / (pointZ-focalZ);
              break;
            case FOC_PRE_CONTINUOUS:
              // Get precomputed focaloffsets in continuous sequential order
              outputPoint[0] = (pointZ*focalX - transformedInputPoint[0]*focalZ
                                + pointZ*(*randomIt)) / (pointZ-focalZ);

              ++randomIt;
              if(randomIt == m_FocalRandomOffsets.end())
                {
                // reset iterator if at end of container
                randomIt = m_FocalRandomOffsets.begin( );
                }

              outputPoint[1] = (pointZ*focalY - transformedInputPoint[1]*focalZ
                                + pointZ*(*randomIt)) / (pointZ-focalZ);

              ++randomIt;
              if(randomIt == m_FocalRandomOffsets.end())
                {
                // reset iterator if at end of container
                randomIt = m_FocalRandomOffsets.begin( );
                }
              break;
            case FOC_PRE_RANDOM:
              // Get precomputed focaloffsets in random order
              randomIt = m_FocalRandomOffsets.begin()
                                + rand.lrand32(m_FocalRandomOffsets.size()-2);
              outputPoint[0] = (pointZ*focalX - transformedInputPoint[0]*focalZ
                                + pointZ*(*randomIt)) / (pointZ-focalZ);

              ++randomIt;
              outputPoint[1] = (pointZ*focalY - transformedInputPoint[1]*focalZ
                                + pointZ*(*randomIt)) / (pointZ-focalZ);
              break;
            case FOC_RANDOM:
              // Compute focaloffsets
              outputPoint[0] = (pointZ*focalX - transformedInputPoint[0]*focalZ
                                + pointZ*rand.normal()*m_FocalRandomOffsetSigma)
                                / (pointZ-focalZ);
              outputPoint[1] = (pointZ*focalY - transformedInputPoint[1]*focalZ
                                + pointZ*rand.normal()*m_FocalRandomOffsetSigma)
                                / (pointZ-focalZ);
              break;
            }

          // Compute corresponding output pixel position
          const bool insideOutputImage = outputPtr->TransformPhysicalPointToContinuousIndex(
              outputPoint, outputContinuousIndex);

          if (insideOutputImage)
            {
            //Transform Pixelvalue
            double inputPixelValueTransformed;
            switch (intensityInterpolation)
              {
              case Superclass::II_NONE:
              default:
                inputPixelValueTransformed = inputPixelValue;
                break;
              case Superclass::II_STEPWISE:
                inputPixelValueTransformed = StepwiseInterpolatePixelValue(
                    intensityTransferFunctionNodes, inputPixelValue);
                break;
              case Superclass::II_LINEAR:
                inputPixelValueTransformed = LinearInterpolatePixelValue(
                    intensityTransferFunctionNodes, inputPixelValue);
                break;
              }

            // Add pixel value to internal outputimage
            switch (m_SplatAssignment)
              {
              case SA_NEAREST:
              default:
                // Convert continuous index to nearest index (nearest neighbour)
                outputIndex.CopyWithRound(outputContinuousIndex);
                // Assign intensity to nearest neighbour
                this->AddToThreadInternalOutputImage(outputIndex,
                    inputPixelValueTransformed, internalIt,
                    outputRegionForThread, threadId);
                break;
              case SA_LINEAR:
                // Assign intensity to the four neighbour pixels
                this->SplatAssignmentLinear(outputContinuousIndex,
                    inputPixelValueTransformed, internalIt,
                    outputRegionForThread, threadId);
                break;
              }
            } // inside image
          } // inside triangle2
        } // inside triangle1
      } // above threshold
    ++inputIt;
    } // while

  return;
}

} // end namespace itk

#endif
