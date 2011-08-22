
#ifndef __ITKPERSPECTIVEPROJECTIONIMAGEFILTER_TXX
#define __ITKPERSPECTIVEPROJECTIONIMAGEFILTER_TXX

#include "itkPerspectiveProjectionImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::PerspectiveProjectionImageFilter()
{
  m_Transform = MatrixOffsetTransformBase<TScalarType,
      itkGetStaticConstMacro(InputImageDimension),
      itkGetStaticConstMacro(InputImageDimension)>::New();
  m_Transform->SetIdentity();
  m_OutputSize.Fill(0);
  m_OutputOrigin.Fill(0.0);
  m_OutputSpacing.Fill(1.0);
  m_FocalPoint.Fill(0.0);
  m_Threshold = 0;
  // m_IntensityTransferFunctionNodes
  m_IntensityInterpolation = II_NONE;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Transform: " << m_Transform.GetPointer() << std::endl;
  os << indent << "OutputSize: " << m_OutputSize << std::endl;
  os << indent << "OutputOrigin: " << m_OutputOrigin << std::endl;
  os << indent << "OutputSpacing: " << m_OutputSpacing << std::endl;
  os << indent << "FocalPoint: " << m_FocalPoint << std::endl;
  os << indent << "Threshold: "
     << static_cast<typename NumericTraits<InputPixelType>::PrintType>(m_Threshold)
     << std::endl;
  os << indent << "IntensityTransferFunctionNodes: ";
  typedef typename IntensityTransferFunctionNodesType::const_iterator ITFIteratorType;
  ITFIteratorType nodeIt = m_IntensityTransferFunctionNodes.begin();
  while (nodeIt != m_IntensityTransferFunctionNodes.end())
    {
    os << "[" << nodeIt->first << ", " << nodeIt->second << "] ";
    ++nodeIt;
    }
  os << std::endl;
  os << indent<< "IntensityInterpolation: " << m_IntensityInterpolation << std::endl;
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::SetOutputSpacing(const double *spacing)
{
  OutputSpacingType s(spacing);
  this->SetOutputSpacing(s);
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::SetOutputOrigin (const OutputOriginType outputOrigin)
{
  itkDebugMacro("setting OutputOrigin to " << outputOrigin);
  if (this->m_OutputOrigin != outputOrigin)
    {
    if(outputOrigin[2] != 0)
      {
      itkExceptionMacro(<< "OutputImage is not in x/y plane (OutputOrigin[2]!=0). " << outputOrigin);
      }

    this->m_OutputOrigin = outputOrigin;
    this->Modified();
    }
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::SetOutputOrigin(const double *origin)
{
  OutputOriginType p(origin);
  this->SetOutputOrigin(p);
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::SetIntensityTransferFunctionNodes(const IntensityTransferFunctionNodesType nodes)
{
  itkDebugMacro("setting IntensityTransferFunctionNodes");
  if (this->m_IntensityTransferFunctionNodes != nodes)
    {
    if(nodes.size() < 2)
      {
      itkExceptionMacro(<< "IntensityTransferFunctionNodes has not enough values to interpolate (<2). " << nodes.size());
      }

    this->m_IntensityTransferFunctionNodes = nodes;
    this->Modified();
    }
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
const typename PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
  ::IntensityTransferFunctionNodesType &
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::GetIntensityTransferFunctionNodes(void) const
{
  itkDebugMacro("returning IntensityTransferFunctionNodes");
  return this->m_IntensityTransferFunctionNodes;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::GenerateOutputInformation(void)
{
  itkDebugMacro("GenerateOutputInformation Start");

  // get pointer to the output
  OutputImagePointer outputPtr = this->GetOutput();
  if (!outputPtr)
    {
    return;
    }

  // Set the size of region
  OutputRegionType outputRegion;
  outputRegion.SetSize(m_OutputSize);
  outputPtr->SetLargestPossibleRegion(outputRegion);

  // Set origin and spacing
  outputPtr->SetOrigin(m_OutputOrigin);
  outputPtr->SetSpacing(m_OutputSpacing);

  itkDebugMacro("GenerateOutputInformation End");
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::GenerateInputRequestedRegion(void)
{
  itkDebugMacro("GenerateInputRequestedRegion Start");

  // call the superclass's implementation of this method
  Superclass::GenerateInputRequestedRegion();

  // get pointer to the input
  if (!this->GetInput())
    {
    return;
    }
  InputImagePointer inputPtr = const_cast<InputImageType*>(this->GetInput());

  // Request the entire input image
  inputPtr->SetRequestedRegionToLargestPossibleRegion();

  itkDebugMacro("GenerateInputRequestedRegion End");
  return;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
unsigned long
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::GetMTime(void) const
{
  unsigned long latestTime = Object::GetMTime();

  if(m_Transform)
    {
    if(latestTime < m_Transform->GetMTime())
      {
      latestTime = m_Transform->GetMTime();
      }
    }

  return latestTime;
}

template <class TInputImage, class TOutputImage, class TScalarType, class TCollector>
void
PerspectiveProjectionImageFilter<TInputImage,TOutputImage,TScalarType,TCollector>
::EnlargeOutputRequestedRegion(DataObject *data)
{
  Superclass::EnlargeOutputRequestedRegion(data);
  data->SetRequestedRegionToLargestPossibleRegion();
  return;
}

} // end namespace itk

#endif
