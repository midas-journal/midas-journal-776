

#ifndef ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_TXX_
#define ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_TXX_

#include <itkMultiResolutionImage2D3DRegistrationMethodCommand.h>
#include <itkImageFileWriter.h>

// support for concrete casting:
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkExhaustiveOptimizerComplete.hxx>


namespace itk
{


template <typename TFramework>
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::MultiResolutionImage2D3DRegistrationMethodCommand()
{
  m_Verbose = false;
  m_Logging = false;
  m_LogFileName = "";
  m_ImageOutput = false;
  m_ImageFileNames.clear();
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::Execute(const Object *, const EventObject &)
{ 
  return;
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::Execute(Object *object, const EventObject &event)
{
  if (!object)
    return;

  // accept the events:
  if(typeid(event) == typeid(StartEvent))
  { 
    this->ExecuteStartEvent(object);
  }
  else if(typeid(event) == typeid(AbortEvent))
  { 
    this->ExecuteAbortEvent(object);
  }
  else if(typeid(event) == typeid(EndEvent))
  { 
    this->ExecuteEndEvent(object);
  }

  return;
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::ExecuteStartEvent(Object *object)
{
  FrameworkPointer fw = SmartPointer<FrameworkType>(
    dynamic_cast<FrameworkType *>(object));

  if (m_Verbose)
    std::cout << "\nThe registration framework has started!" << std::endl;

  m_StartTimeStamp = fw->GetClock()->GetTimeStamp(); // store start time

  return;
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::ExecuteAbortEvent(Object *object)
{
  if (m_Verbose)
    std::cout << "\nThe registration framework has been aborted!" <<
      std::endl;

  return;
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::ExecuteEndEvent(Object *object)
{
  if (m_Verbose)
    std::cout << "\nThe registration framework has finished!" << std::endl;

  FrameworkPointer fw = SmartPointer<FrameworkType>(
    dynamic_cast<FrameworkType *>(object));

  // get number of optimization iterations, final transformation parameters
  // and best metric value:
  unsigned int numberOfIterations;
  typename FrameworkType::BaseOptimizerType::ParametersType finalPars;
  double bestMetric;

  // get current registration results:
  if (GetCurrentResults(fw, finalPars, bestMetric, numberOfIterations))
  {
    // verbose and log final registration results:
    VerboseAndLogFinalRegistrationResult(fw, numberOfIterations, finalPars,
      bestMetric);

    // write out final moving and PDF images:
    WriteFinalMovingImage(fw, finalPars);
    WriteFinalPDFImage(fw, finalPars);
    
    // be sure that the final transform parameters are set:
    fw->GetTransformation()->SetParameters(finalPars);
  }

  return;
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::WriteFinalMovingImage(FrameworkPointer fw,
  typename FrameworkType::BaseOptimizerType::ParametersType finalPars)
{
  if (m_ImageOutput && (m_ImageFileNames.size() >= 1) &&
      (m_ImageFileNames[0].length() > 0))
  { 
    if (m_Verbose)
      std::cout << std::endl << "Writing final moving image to disc." <<
        std::endl << std::endl;

    // to get the final moving image for current level, the transformation
    // must be set to the result transformation parameters:
    typedef typename FrameworkType::BaseTransformType::ParametersType
      ParametersType;
    ParametersType parsStore = fw->GetTransformation()->GetParameters();

    fw->GetTransformation()->SetParameters(finalPars);

    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    // - final moving image:
    char fn[2048];
    sprintf(fn, m_ImageFileNames[0].c_str(),
      fw->GetRegistration()->GetNumberOfLevels());
    writer->SetFileName(fn);

    typename ImageType::ConstPointer mov = fw->GetCurrentProjectionImage();
    if (mov)
    { 
      writer->SetInput(mov);
      try
      { 
        writer->Update();
      }
      catch (ExceptionObject &e)
      { 
        std::cerr << "Error during writing final moving image (" <<
          fn << "): " << e << std::endl;
      }
    }

    // restore parameters:
    fw->GetTransformation()->SetParameters(parsStore);
  }
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::WriteFinalPDFImage(FrameworkPointer fw,
  typename FrameworkType::BaseOptimizerType::ParametersType finalPars)
{
  if (m_ImageOutput && (m_ImageFileNames.size() >= 2) &&
      (m_ImageFileNames[1].length() > 0))
  { 
    if (m_Verbose)
      std::cout << std::endl << "Writing final PDF image to disc." <<
        std::endl << std::endl;

    // to get the final PDF image for current level, the transformation
    // must be set to the result transformation parameters:
    typedef typename FrameworkType::BaseTransformType::ParametersType
      ParametersType;
    ParametersType parsStore = fw->GetTransformation()->GetParameters();

    fw->GetTransformation()->SetParameters(finalPars);

    typedef typename FrameworkType::HistImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    // - final PDF image:
    char fn[2048];
    sprintf(fn, m_ImageFileNames[1].c_str(),
      fw->GetRegistration()->GetNumberOfLevels());
    writer->SetFileName(fn);

    typename ImageType::Pointer pdf = fw->GetCurrentMetricHistogram(true);

    if (pdf)
    { 
      writer->SetInput(pdf);
      try
      { 
        writer->Update();
      }
      catch (ExceptionObject &e)
      { 
        std::cerr << "Error during writing final PDF image (" <<
          fn << "): " << e << std::endl;
      }
    }

    // restore parameters:
    fw->GetTransformation()->SetParameters(parsStore);
  }
}

template <typename TFramework>
void 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::VerboseAndLogFinalRegistrationResult(FrameworkPointer fw,
  unsigned int numberOfIterations,
  typename FrameworkType::BaseOptimizerType::ParametersType finalPars,
  double bestValue)
{
  std::vector<std::string> vs;
  std::ostringstream os;

  vs.push_back("\n\nInital Transformation:\n");
  os.str(""); os << " (*) Paramters: \t\t" << 
    fw->GetInitialTransformationParamters();
  vs.push_back(os.str());
  
  vs.push_back("\n\nFinal Registration Parameters:\n");
  os.str(""); os << " (*) Number of Iterations: \t\t" << numberOfIterations;
  vs.push_back(os.str());
  os.str(""); os << " (*) Final Parameters: \t\t" << finalPars;
  vs.push_back(os.str());
  os.str(""); os << " (*) Best Metric: \t\t" << bestValue;
  vs.push_back(os.str());
  typename FrameworkType::RealTimeStampType totalTime =
    fw->GetClock()->GetTimeStamp() - m_StartTimeStamp;
  os.str(""); os << " (*) Total Time: \t\t" << totalTime << "s";
  vs.push_back(os.str());

  if (m_Verbose)
    for (unsigned int i = 0; i < vs.size(); ++i)
      std::cout << vs[i] << std::endl;

  if (m_Logging && (m_LogFileName.length() > 0))
  { 
    char fn[2048];

    sprintf(fn, m_LogFileName.c_str(),
      fw->GetRegistration()->GetNumberOfLevels());
    m_LogFile.open(fn);

    for (unsigned int i = 0; i < vs.size(); ++i)
      m_LogFile << vs[i] << std::endl;

    m_LogFile.close();
  }
}

template <typename TFramework>
bool 
MultiResolutionImage2D3DRegistrationMethodCommand<TFramework>
::GetCurrentResults(typename FrameworkType::Pointer fw,
  typename FrameworkType::BaseOptimizerType::ParametersType &finalPars,
  double &bestValue, unsigned int &numberOfIterations)
{
  bool validResults = false;

  typedef ExhaustiveOptimizerComplete ExhOptimizerType;
  typedef ExhOptimizerType::ConstPointer ExhOptimizerConstPointer;
  typedef OnePlusOneEvolutionaryOptimizer EvolOptimizerType;
  typedef EvolOptimizerType::ConstPointer EvolOptimizerConstPointer;

  // do concrete casts (first these are pure dummies, then casted):
  ExhOptimizerConstPointer exhOptimizer = SmartPointer<const ExhOptimizerType>(
    dynamic_cast<const ExhOptimizerType *>(
    ExhOptimizerType::New().GetPointer()));
  EvolOptimizerConstPointer evolOptimizer = 
    SmartPointer<const EvolOptimizerType>(
    dynamic_cast<const EvolOptimizerType *>(
    EvolOptimizerType::New().GetPointer()));

  if (fw->GetOptimizer()->GetNameOfClass() == exhOptimizer->GetNameOfClass())
  { 
    exhOptimizer =
      SmartPointer<const ExhOptimizerType>(
      dynamic_cast<const ExhOptimizerType *>(
      fw->GetOptimizer().GetPointer()));

    numberOfIterations = exhOptimizer->GetCurrentIteration() + 1;

    if (fw->GetOptimizationMinimization())
    { 
      finalPars = exhOptimizer->GetMinimumMetricValuePosition();
      bestValue = exhOptimizer->GetMinimumMetricValue();
    }
    else
    { 
      finalPars = exhOptimizer->GetMaximumMetricValuePosition();
      bestValue = exhOptimizer->GetMaximumMetricValue();
    }
    validResults = true;
  }
  else if (fw->GetOptimizer()->GetNameOfClass() == 
           evolOptimizer->GetNameOfClass())
  { 
    evolOptimizer =
      SmartPointer<const EvolOptimizerType>(
      dynamic_cast<const EvolOptimizerType *>(
      fw->GetOptimizer().GetPointer()));

    numberOfIterations = evolOptimizer->GetCurrentIteration();
    finalPars = evolOptimizer->GetCurrentPosition();
    bestValue = evolOptimizer->GetCurrentCost();

    validResults = true;
  }

  return validResults;
}


}


#endif /* ITKMULTIRESOLUTIONIMAGE2D3DREGISTRATIONMETHODCOMMAND_TXX_ */
