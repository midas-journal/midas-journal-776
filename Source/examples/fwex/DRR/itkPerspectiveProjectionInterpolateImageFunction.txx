
#ifndef __ITKPERSPECTIVEPROJECTIONINTERPOLATEIMAGEFUNCTION_TXX
#define __ITKPERSPECTIVEPROJECTIONINTERPOLATEIMAGEFUNCTION_TXX

#include "itkPerspectiveProjectionInterpolateImageFunction.h"
#include "itkMinimumMaximumImageCalculator.h"

namespace itk
{

template<class TInputImage, class TCoordRep, class TCollector>
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::PerspectiveProjectionInterpolateImageFunction()
{
  m_PerspectiveProjectionImageFilter = NULL;
  m_ProjectedImage = NULL;
  m_ProjectionInterpolator = NearestNeighborInterpolateImageFunction<ProjectionImageType, TCoordRep>::New();;
  m_InverseTransform = InverseTransformType::New();
  m_UpdateCommand = UpdateCommandType::New();
  m_UpdateCommand->SetTarget(this);
  m_UpdateCommand->SetUpdateActive(true);
  m_UpdateCommandTag = 0;

  m_FindMinMaxIntensityMode = false; // auto-search min/max intensities
  m_MinIntensity = NumericTraits<
    typename PerspectiveProjectionImageFilterType::OutputPixelType>::max();
  m_MaxIntensity = NumericTraits<
    typename PerspectiveProjectionImageFilterType::OutputPixelType>::min();
  m_IntermediateFilter = NULL; // no filter is applied
}

template<class TInputImage, class TCoordRep, class TCollector>
void
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);

  os << indent << "PerspectiveProjectionImageFilter: "
     << m_PerspectiveProjectionImageFilter.GetPointer() << std::endl;
  os << indent << "ProjectedImage: " << m_ProjectedImage.GetPointer() << std::endl;
  os << indent << "InternalInterpolator: " << m_ProjectionInterpolator.GetPointer() << std::endl;
  os << indent << "InverseTransform: " << m_InverseTransform.GetPointer() << std::endl;
  os << indent << "UpdateCommand: " << m_UpdateCommand.GetPointer() << std::endl;
  os << indent<< "UpdateCommandTag: "
     << static_cast<typename NumericTraits<unsigned long>::PrintType>(m_UpdateCommandTag)
     << std::endl;
  return;
}

template<class TInputImage, class TCoordRep, class TCollector>
void
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::Update()
{
  // update projected image
  this->m_PerspectiveProjectionImageFilter->Update();
  // basically the raw projected image is set (BUT: an optional intermediate
  // filter can be applied)
  this->m_ProjectedImage = this->m_PerspectiveProjectionImageFilter->GetOutput();

  // apply the optional intermediate filter if defined
  if (this->m_IntermediateFilter)
    {
    this->m_IntermediateFilter->SetInput(
      this->m_PerspectiveProjectionImageFilter->GetOutput());
    try
      {
      // (the filter is expected to be externally pre-configured)
      this->m_IntermediateFilter->UpdateLargestPossibleRegion();
      this->m_IntermediateFilter->Update();
      this->m_ProjectedImage = this->m_IntermediateFilter->GetOutput();
      }
    catch (ExceptionObject &e)
      {
      itkExceptionMacro(<< "ERROR (during intermediate filtering): " << e);
      }
    }

  // find min/max intensity values if demanded:
  if (m_FindMinMaxIntensityMode)
    {
    typedef MinimumMaximumImageCalculator<ProjectionImageType> MinMaxCalcType;

    typename MinMaxCalcType::Pointer minMaxCalc = MinMaxCalcType::New();
    minMaxCalc->SetImage(this->m_ProjectedImage);

    minMaxCalc->Compute();

    m_MinIntensity = minMaxCalc->GetMinimum();
    m_MaxIntensity = minMaxCalc->GetMaximum();
    }

  this->m_ProjectionInterpolator->SetInputImage(m_ProjectedImage);

  // update inverse transform
  InverseTransformTypePointer transform = dynamic_cast<InverseTransformType *>(
      this->m_PerspectiveProjectionImageFilter->GetTransform());
  m_InverseTransform->SetCenter(transform->GetCenter());
  transform->GetInverse(m_InverseTransform);

  this->Modified();
  return;
}

template<class TInputImage, class TCoordRep, class TCollector>
void
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::SetPerspectiveProjectionImageFilter (PerspectiveProjectionImageFilterType* perspectiveProjectionImageFilter)
{
  itkDebugMacro("setting PerspectiveProjectionImageFilter to " << perspectiveProjectionImageFilter);
  if (this->m_PerspectiveProjectionImageFilter != perspectiveProjectionImageFilter)
    {
    if (this->m_PerspectiveProjectionImageFilter)
      {
      // remove observer first
      this->m_PerspectiveProjectionImageFilter->GetTransform()->RemoveObserver(
          this->m_UpdateCommandTag);
      this->m_UpdateCommandTag = 0;
      }

    // set new projection
    this->m_PerspectiveProjectionImageFilter = perspectiveProjectionImageFilter;

    // add observer
    if (this->m_PerspectiveProjectionImageFilter &&
        this->m_PerspectiveProjectionImageFilter->GetTransform())
      {
      this->m_UpdateCommandTag = this->m_PerspectiveProjectionImageFilter->GetTransform()
        ->AddObserver(ModifiedEvent(), this->m_UpdateCommand);
      }

    // force update
    this->Update();
    }
  return;
}

template<class TInputImage, class TCoordRep, class TCollector>
void
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::SetUpdateActive (bool updateActive)
{
  itkDebugMacro("setting UpdateActive to " << updateActive);
  if (this->m_UpdateCommand->GetUpdateActive() != updateActive)
    {
    m_UpdateCommand->SetUpdateActive(updateActive);
    this->Update();
    }
  return;
}

template<class TInputImage, class TCoordRep, class TCollector>
bool
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::GetUpdateActive (void) const
{
  return this->m_UpdateCommand->GetUpdateActive();
}

template<class TInputImage, class TCoordRep, class TCollector>
void
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::SetInputImage(const InputImageType * ptr)
{
  Superclass::SetInputImage(ptr);
  this->m_PerspectiveProjectionImageFilter->SetInput(ptr);
  this->Modified();
  return;
}

template<class TInputImage, class TCoordRep, class TCollector>
typename PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep,
  TCollector>::OutputType
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::Evaluate(const PointType& point) const
{
  PointType transformedPoint = m_InverseTransform->TransformPoint(point);
  // no conversion is needed: point is in physical space
  // convert point to index in projected image and get intensity
  ContinuousIndexType cindex;
  this->m_ProjectedImage->TransformPhysicalPointToContinuousIndex(transformedPoint, cindex);
  return (this->m_ProjectionInterpolator->EvaluateAtContinuousIndex(cindex));
}

template<class TInputImage, class TCoordRep, class TCollector>
typename PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep,
  TCollector>::OutputType
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::EvaluateAtContinuousIndex(const ContinuousIndexType & index) const
{
  // conversion is needed: index is in inputimage (moving image) and not in DRR
  // convert index to point
  PointType point;
  this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(index, point);
  return (this->Evaluate(point));
}

template<class TInputImage, class TCoordRep, class TCollector>
typename PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep,
  TCollector>::OutputType
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::EvaluateAtIndex(const IndexType & index) const
{
  // conversion is needed: index is in inputimage (moving image) and not in DRR
  // convert index to point
  PointType point;
  this->GetInputImage()->TransformIndexToPhysicalPoint(index, point);
  return (this->Evaluate(point));
}

template<class TInputImage, class TCoordRep, class TCollector>
bool
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::IsInsideBuffer(const PointType & point) const
{
  PointType transformedPoint = m_InverseTransform->TransformPoint(point);
  // no conversion is needed: point is in physical space
  // convert point to index in projected image and check if inside buffer
  ContinuousIndexType cindex;
  this->m_ProjectedImage->TransformPhysicalPointToContinuousIndex(transformedPoint, cindex);
  IndexType nindex;
  this->ConvertContinuousIndexToNearestIndex(cindex, nindex);
  return (this->m_ProjectionInterpolator->IsInsideBuffer(nindex));
}

template<class TInputImage, class TCoordRep, class TCollector>
bool
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::IsInsideBuffer(const IndexType & index) const
{
  // conversion is needed: index is in inputimage (moving image) and not in DRR
  // convert index to point
  PointType point;
  this->GetInputImage()->TransformIndexToPhysicalPoint(index, point);
  return (this->IsInsideBuffer(point));
}

template<class TInputImage, class TCoordRep, class TCollector>
bool
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::IsInsideBuffer(const ContinuousIndexType & index) const
{
  // conversion is needed: index is in inputimage (moving image) and not in DRR
  // convert index to point
  PointType point;
  this->GetInputImage()->TransformContinuousIndexToPhysicalPoint(index, point);
  return (this->IsInsideBuffer(point));
}

template<class TInputImage, class TCoordRep, class TCollector>
unsigned long
PerspectiveProjectionInterpolateImageFunction<TInputImage, TCoordRep, TCollector>
::GetMTime(void) const
{
  unsigned long latestTime = Object::GetMTime();

  if(m_PerspectiveProjectionImageFilter)
    {
    if(latestTime < m_PerspectiveProjectionImageFilter->GetMTime())
      {
      latestTime = m_PerspectiveProjectionImageFilter->GetMTime();
      }
    }

  if(m_ProjectionInterpolator)
    {
    if( latestTime < m_ProjectionInterpolator->GetMTime() )
      {
      latestTime = m_ProjectionInterpolator->GetMTime();
      }
    }

  return latestTime;
}

} // end namespace itk

#endif
