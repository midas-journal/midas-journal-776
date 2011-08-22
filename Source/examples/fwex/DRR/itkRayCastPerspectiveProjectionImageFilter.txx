
#ifndef __ITKRAYCASTPERSPECTIVEPROJECTIONIMAGEFILTER_TXX
#define __ITKRAYCASTPERSPECTIVEPROJECTIONIMAGEFILTER_TXX

#include "itkRayCastPerspectiveProjectionImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkContinuousIndex.h"
#include "itkRealTimeClock.h"
#include "itkPoint.h"
#include "itkVector.h"
#include "itkCovariantVector.h"
#include "itkFixedArray.h"
#include "vnl/vnl_random.h"
#include "vcl_cmath.h"
#include "vnl/vnl_math.h"
#include <algorithm>

namespace itk
{

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::RayCastPerspectiveProjectionImageFilter()
{
  m_InverseTransform = MatrixOffsetTransformBase<TScalarType,
      itkGetStaticConstMacro(InputImageDimension),
      itkGetStaticConstMacro(InputImageDimension)>::New();
  this->GetTransform()->GetInverse(m_InverseTransform);
  m_Interpolator = LinearInterpolateImageFunction<InputImageType, TScalarType>::New();
  m_RayStepSize = 1;
  m_InternalOutputImage = InternalOutputImageType::New();
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "InverseTransform: " << m_InverseTransform.GetPointer() << std::endl;
  os << indent << "Interpolator: " << m_Interpolator.GetPointer() << std::endl;
  os << indent << "RayStepSize: "
     << static_cast<typename NumericTraits<TransformScalarType>::PrintType>(m_RayStepSize)
     << std::endl;
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::BeforeThreadedGenerateData(void)
{
  // Get members from superclass
  TransformTypeConstPointer transform = this->GetTransform();
  const IntensityTransferFunctionNodesType& intensityTransferFunctionNodes =
      this->GetIntensityTransferFunctionNodes();
  const IntensityInterpolationEnumType& intensityInterpolation =
      this->GetIntensityInterpolation();

  // Check if everything is set properly
  if(!transform)
    {
    itkExceptionMacro(<< "Transform not set");
    }

  if(intensityInterpolation != Superclass::II_NONE && intensityTransferFunctionNodes.size() < 2)
    {
    itkExceptionMacro(<< "IntensityTransferFunctionNodes has not enough values to interpolate (<2).");
    }

  if(!m_InverseTransform)
    {
    itkExceptionMacro(<< "InverseTransform not set");
    }

  if(!m_Interpolator)
    {
    itkExceptionMacro(<< "Interpolator not set");
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

  // Get inverse transform
  m_InverseTransform->SetCenter(transform->GetCenter());
  transform->GetInverse(m_InverseTransform);

  // Connect input image to interpolator
  m_Interpolator->SetInputImage(this->GetInput());

  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::AfterThreadedGenerateData(void)
{
  // Disconnect input image from the interpolator
  m_Interpolator->SetInputImage(NULL);

  // Get the output pointer
  OutputImagePointer outputPtr = this->GetOutput();

  // minimum/maximum values of the OutputPixelType
  const OutputPixelType outputMinPixelValue = NumericTraits<OutputPixelType>::NonpositiveMin();
  const OutputPixelType outputMaxPixelValue = NumericTraits<OutputPixelType>::max();
  const CollectorValueType minPixelValue = static_cast<CollectorValueType>(outputMinPixelValue);
  const CollectorValueType maxPixelValue = static_cast<CollectorValueType>(outputMaxPixelValue);

  // Create iterator to access internal outputimage
  typedef ImageRegionIterator<InternalOutputImageType> InternalIteratorType;
  InternalIteratorType internalIt = InternalIteratorType(m_InternalOutputImage,
      m_InternalOutputImage->GetLargestPossibleRegion());

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
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::ThreadedGenerateData(const OutputRegionType &outputRegionForThread, int threadId)
{
  // Get the input and output ConstPointer
  InputImageConstPointer inputPtr = this->GetInput();
  OutputImageConstPointer outputPtr = this->GetOutput();

  // Get members from superclass
  TransformTypeConstPointer transform = this->GetTransform();
  // const OutputSizeType& outputSize = this->GetOutputSize();
  // const OutputOriginType& outputOrigin = this->GetOutputOrigin();
  // const OutputSpacingType& outputSpacing = this->GetOutputSpacing();
  const InputPointType& focalPoint = this->GetFocalPoint();
  const InputPixelType& threshold = this->GetThreshold();
  const IntensityTransferFunctionNodesType& intensityTransferFunctionNodes =
      this->GetIntensityTransferFunctionNodes();
  const IntensityInterpolationEnumType& intensityInterpolation =
      this->GetIntensityInterpolation();

  // Create iterator to access internal outputimage
  typedef ImageRegionIteratorWithIndex<InternalOutputImageType> InternalIteratorType;
  InternalIteratorType internalIt = InternalIteratorType(m_InternalOutputImage,
      outputRegionForThread);

  // Define a few variables that will be used for sampling along the rays

  // constant references to step size
  const double& rayStepSize = m_RayStepSize;

  // Cartesian coordinates of current output pixel (internal image)
  OutputPointType outputPoint;
  outputPoint.Fill(0);

  // Direction of the ray that passes through the focal point and the output pixel
  typedef Vector <double,itkGetStaticConstMacro(InputImageDimension)> VectorType;
  VectorType directionVector;
  directionVector.Fill(0);

  // Cartesian coordinates of current sample point on ray that is moved by step size
  InputPointType rayStepPoint;
  rayStepPoint.Fill(0);

  // Continuous index of current sample point on ray that is moved with step size
  typedef ContinuousIndex<typename InputPointType::ValueType,
      itkGetStaticConstMacro(InputImageDimension)> InputContinuousIndexType;
  InputContinuousIndexType rayStepContinuousIndex;
  rayStepContinuousIndex.Fill(0);

  // Define a few points and vectors to check if a ray hits the volume

  // Convert point to vector to perform plane intersection calculations
  VectorType outputPointVector;
  outputPointVector.SetVnlVector(outputPoint.GetVnlVector());

  // Get indices of corners from volume
  const InputRegionType inputRegion = inputPtr->GetLargestPossibleRegion();
  typedef typename InputRegionType::IndexType InputRegionIndexType;
  typedef typename InputRegionType::SizeType InputRegionSizeType;
  typedef typename InputRegionType::IndexValueType InputRegionIndexValueType;

  // smallest index of input image
  const InputRegionIndexType inputRegionIndexMin = inputRegion.GetIndex();
  const InputRegionSizeType inputRegionSize = inputRegion.GetSize();

  // highest index of input image
  InputRegionIndexType inputRegionIndexMax;
  inputRegionIndexMax[0] = inputRegionIndexMin[0] +
      static_cast<InputRegionIndexValueType>(inputRegionSize[0]) - 1;
  inputRegionIndexMax[1] = inputRegionIndexMin[1] +
      static_cast<InputRegionIndexValueType>(inputRegionSize[1]) - 1;
  inputRegionIndexMax[2] = inputRegionIndexMin[2] +
      static_cast<InputRegionIndexValueType>(inputRegionSize[2]) - 1;

  // Vertices of input image region in input image space and transformed space
  InputPointType inputRegionPointMin;
  inputPtr->TransformIndexToPhysicalPoint(inputRegionIndexMin, inputRegionPointMin);
  const InputPointType inputRegionPointMinTransformed =
      transform->TransformPoint(inputRegionPointMin);

  InputPointType inputRegionPointMax;
  inputPtr->TransformIndexToPhysicalPoint(inputRegionIndexMax, inputRegionPointMax);
  const InputPointType inputRegionPointMaxTransformed =
      transform->TransformPoint(inputRegionPointMax);

  // Convert points to vectors to perform plane intersection calculations
  VectorType inputRegionPointMinTransformedVector;
  inputRegionPointMinTransformedVector.SetVnlVector(inputRegionPointMinTransformed.GetVnlVector());
  VectorType inputRegionPointMaxTransformedVector;
  inputRegionPointMaxTransformedVector.SetVnlVector(inputRegionPointMaxTransformed.GetVnlVector());

  // Normals of planes in input image space and transformed space
  typedef CovariantVector <double,
      itkGetStaticConstMacro(InputImageDimension)> CovariantVectorType;
  CovariantVectorType xPlaneNormalVector;
  xPlaneNormalVector.Fill(0);
  xPlaneNormalVector[1] = 1;
  const CovariantVectorType xPlaneNormalVectorTransformed =
      transform->TransformCovariantVector(xPlaneNormalVector);

  CovariantVectorType yPlaneNormalVector;
  yPlaneNormalVector.Fill(0);
  yPlaneNormalVector[0] = 1;
  const CovariantVectorType yPlaneNormalVectorTransformed =
      transform->TransformCovariantVector(yPlaneNormalVector);

  CovariantVectorType zPlaneNormalVector;
  zPlaneNormalVector.Fill(0);
  zPlaneNormalVector[2] = 1;
  const CovariantVectorType zPlaneNormalVectorTransformed =
      transform->TransformCovariantVector(zPlaneNormalVector);

  // Container for computed intersections of ray with planes at k
  typedef FixedArray<double,6> FixedArrayType;
  FixedArrayType kIntersections;
  kIntersections.Fill(0);

  // Container for computed intersections of ray with volume at k
  FixedArrayType kVolumeIntersections;
  kVolumeIntersections.Fill(0);

  // Cartesian coordinates of volume intersection
  InputPointType intersectionPoint;
  intersectionPoint.Fill(0);

  // Cartesian coordinates of volume intersection in input image space
  InputPointType intersectionPointTransformed;
  intersectionPointTransformed.Fill(0);

  // Index of volume intersection in input image
  InputIndexType intersectionIndex;
  intersectionIndex.Fill(0);

  InputContinuousIndexType intersectionContinuousIndex;
  intersectionContinuousIndex.Fill(0);

  // Iterate over output image
  internalIt.GoToBegin();

  while (!internalIt.IsAtEnd())
    {
    // Get the coordinates of the output pixel
    outputPtr->TransformIndexToPhysicalPoint(internalIt.GetIndex(), outputPoint);

    // direction vector that passes through pixel and focalpoint
    directionVector = focalPoint - outputPoint;

    // steps needed between pixel and focalpoint
    typedef int KValueType;
    const KValueType kMax = static_cast<KValueType>(directionVector.GetNorm()/rayStepSize);

    // normalise and scale direction vector to step size
    directionVector.Normalize();
    directionVector *= rayStepSize;

    // Find plane intersections to limit ray traversal to volume
    // calculate intersections of ray with volume
    // 6 possible intersections because of 6 planes
    outputPointVector.SetVnlVector(outputPoint.GetVnlVector());

    // intersection of ray with x-planes (min- and max-plane)
    double denominator = xPlaneNormalVectorTransformed * directionVector;
    if(denominator != 0)
      {
      kIntersections[0] = (xPlaneNormalVectorTransformed * inputRegionPointMinTransformedVector
          - xPlaneNormalVectorTransformed * outputPointVector) / denominator;
      kIntersections[1] = (xPlaneNormalVectorTransformed * inputRegionPointMaxTransformedVector
          - xPlaneNormalVectorTransformed * outputPointVector) / denominator;
      }
    else
      {
      // if division by zero
      kIntersections[0] = NumericTraits<KValueType>::max();
      kIntersections[1] = NumericTraits<KValueType>::max();
      }

    // intersection of ray with y-planes (min and max-plane)
    denominator = yPlaneNormalVectorTransformed * directionVector;
    if(denominator != 0)
      {
      kIntersections[2] = (yPlaneNormalVectorTransformed * inputRegionPointMinTransformedVector
          - yPlaneNormalVectorTransformed * outputPointVector) / denominator;
      kIntersections[3] = (yPlaneNormalVectorTransformed * inputRegionPointMaxTransformedVector
          - yPlaneNormalVectorTransformed * outputPointVector)/ denominator;
      }
    else
      {
      // if division by zero
      kIntersections[2] = NumericTraits<KValueType>::max();
      kIntersections[3] = NumericTraits<KValueType>::max();
      }

    // intersection of ray with z-planes (min and max-plane)
    denominator = zPlaneNormalVectorTransformed * directionVector;
    if(denominator != 0)
      {
      kIntersections[4] = (zPlaneNormalVectorTransformed * inputRegionPointMinTransformedVector
          - zPlaneNormalVectorTransformed * outputPointVector) / denominator;
      kIntersections[5] = (zPlaneNormalVectorTransformed * inputRegionPointMaxTransformedVector
          - zPlaneNormalVectorTransformed * outputPointVector) / denominator;
      }
    else
      {
      // if division by zero
      kIntersections[4] = NumericTraits<KValueType>::max();
      kIntersections[5] = NumericTraits<KValueType>::max();
      }

    // check if a ray intersects with volume
    kVolumeIntersections.Fill(0);
    int intersectionCounter = 0;

    // check intersections if they hit volume
    for (FixedArrayType::ConstIterator it = kIntersections.Begin(); it != kIntersections.End(); ++it)
      {
      const double& k = *it;
      // compute coordinates of intersection point
      intersectionPoint = outputPoint + directionVector * k;
      // transform point to input image space
      intersectionPointTransformed = m_InverseTransform->TransformPoint(intersectionPoint);
      inputPtr->TransformPhysicalPointToContinuousIndex(
          intersectionPointTransformed,intersectionContinuousIndex);
      // rounding to nearest index because of rounding errors
      // Voronoi region is used to check if intersection point hits volume
      intersectionIndex.CopyWithRound(intersectionContinuousIndex);
      if (m_Interpolator->IsInsideBuffer(intersectionIndex))
        {
        // if intersection point hits volume it is saved
        if(k > NumericTraits<KValueType>::max())
          {
          kVolumeIntersections[intersectionCounter] = NumericTraits<KValueType>::max();
          }
        else
          {
          kVolumeIntersections[intersectionCounter] = k;
          }
        intersectionCounter++;
        }
      }

    // check if ray hits volume
    if (intersectionCounter >= 2)
      {
      // Ray hits volume if it has two or more (Voronoi region is used because of
      // rounding errors) intersection points
      // Sort intersection point positions on ray and choose smallest and highest
      std::sort(kVolumeIntersections.Begin(), kVolumeIntersections.End());
      KValueType kFirst = static_cast<KValueType>(vcl_floor(
          kVolumeIntersections[kVolumeIntersections.Length - intersectionCounter]));
      KValueType kLast = static_cast<KValueType>(vcl_ceil(
          kVolumeIntersections[kVolumeIntersections.Length - 1]));
      kFirst = vnl_math_max(NumericTraits<KValueType>::ZeroValue(), kFirst);
      kLast = vnl_math_min(kMax, kLast);

      // Traverse ray
      for (KValueType k = kFirst; k <= kLast; ++k)
        {
        // Compute corresponding ray step position
        rayStepPoint = outputPoint + (directionVector * k);

        // Transform point to input image space and convert it to continuous index
        rayStepPoint = m_InverseTransform->TransformPoint(rayStepPoint);
        const bool insideInputImage = inputPtr->TransformPhysicalPointToContinuousIndex(
            rayStepPoint, rayStepContinuousIndex);

        // Check if point is inside the input image
        if (insideInputImage)
          {
          // Interpolate intensity
          const double rayStepPixelValue = m_Interpolator->EvaluateAtContinuousIndex(
              rayStepContinuousIndex);

          // Check if pixelvalue is above threshold
          if (rayStepPixelValue >= threshold)
            {
            //Transform Pixelvalue
            double rayStepPixelValueTransformed;
            switch (intensityInterpolation)
              {
              case Superclass::II_NONE:
              default:
                rayStepPixelValueTransformed = rayStepPixelValue;
                break;
              case Superclass::II_STEPWISE:
                rayStepPixelValueTransformed = StepwiseInterpolatePixelValue(
                    intensityTransferFunctionNodes, rayStepPixelValue);
                break;
              case Superclass::II_LINEAR:
                rayStepPixelValueTransformed = LinearInterpolatePixelValue(
                    intensityTransferFunctionNodes, rayStepPixelValue);
                break;
              }

            // Add pixel value to internal outputimage
            internalIt.Value().Add(rayStepPixelValueTransformed);
            } // above threshold
          } // inside image
        } // ray traversal
      }
    else
      {
      // intersectionCounter == 0: Ray does not hit the volume
      // intersectionCounter == 1: Ray hits exactly on an edge or corner
      // due to rounding errors this case is not handled instead
      // the Voronoi region is used to avoid rounding errors
      } // volume intersection
    ++internalIt;
    } // while

  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
unsigned long
RayCastPerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::GetMTime(void) const
{
  unsigned long latestTime = Object::GetMTime();

  TransformTypeConstPointer transform = this->GetTransform();

  if(transform)
    {
    if(latestTime < transform->GetMTime())
      {
      latestTime = transform->GetMTime();
      }
    }

  if( m_Interpolator )
    {
    if( latestTime < m_Interpolator->GetMTime() )
      {
      latestTime = m_Interpolator->GetMTime();
      }
    }

  return latestTime;
}

} // end namespace itk

#endif
