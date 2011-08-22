

#ifndef ITKOPTIMIZERITERATIONCOMMAND_TXX_
#define ITKOPTIMIZERITERATIONCOMMAND_TXX_

#include <itkOptimizerIterationCommand.h>

#include <itkImageFileWriter.h>
#include <itkMath.h>

// support for concrete casting:
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkExhaustiveOptimizerComplete.hxx>

namespace itk
{


template <typename TFramework>
OptimizerIterationCommand<TFramework>
::OptimizerIterationCommand()
{
  m_Clock = itk::RealTimeClock::New();
  m_LogFileName = "";
  m_Logging = false;
  m_LogCurrentParametersAlso = false;
  m_LogSplitValue = -1;
  m_CurrLogSplitId = 0;
  m_ImageOutput = false;
  m_ImageAutoOutput = false;
  m_ImageModulo = 10;
  m_ImageBaseFileName = "";
  m_PDFOutput = false;
  m_PDFAutoOutput = false;
  m_PDFModulo = 10;
  m_PDFBaseFileName = "";
  m_Framework = NULL;
  m_LastCostValue = -1; // it's more or less uninteresting ...
}

template <typename TFramework>
OptimizerIterationCommand<TFramework>
::~OptimizerIterationCommand()
{ 
  if ((m_LogFileName.length() > 0) && m_Logging)
    m_LogFile.close();
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::Execute(Object *caller, const EventObject &event)
{
  Execute((const itk::Object *)caller, event);
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::Execute(const Object *object, const EventObject &event)
{
  if (!m_Framework || !object)
    return;

  // accept iterations only!
  if(typeid(event) != typeid(IterationEvent))
  { 
    return;
  }

  int currIt = -1;
  int maxIt = -1;
  typename BaseOptimizerType::ParametersType pos;
  typename BaseMetricType::MeasureType thisIterationCostValue;

  // cast the basic components:
  BaseOptimizerConstPointer optimizer = SmartPointer<const BaseOptimizerType>(
    dynamic_cast<const BaseOptimizerType *>(object));

  // get the parameters of current iteration:
  if (GetCurrentIterationParameters(optimizer, currIt, maxIt,
      thisIterationCostValue, pos))
  {
    // open/close/split central log-file:
    OpenCloseSplitLogFile(currIt);

    // verbose and log the registration parameters at current iteration:
    VerboseAndLogRegistrationParameters(currIt, thisIterationCostValue, pos);

    // verbose registration progress information:
    VerboseProgress(currIt, maxIt);


    // check whether cost value changed since last iteration:
    bool costValueChanged = ((currIt == 1) ||
      ((currIt > 1) && (thisIterationCostValue != m_LastCostValue)));
    if (costValueChanged)
      m_LastCostValue = thisIterationCostValue; // store!

    // write out current moving and PDF images:
    WriteMovingImage(currIt, costValueChanged);
    WritePDFImage(currIt, costValueChanged);
  }
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::VerboseProgress(int currIt, int maxIt)
{
  if (m_Verbose)
  { 
    float prog;
    std::ostringstream leftTimeStr;

    prog = (float)currIt / (float)maxIt * 100;

    leftTimeStr.str("");
    leftTimeStr << "<unknown>";

    if (currIt == 1)
      m_StartTime = m_Clock->GetTimeStamp();

    // progress/left time estimation:
    if (((currIt % 20) == 0) || (maxIt == currIt))
    { 
      if (currIt > 1)
      { 
        double now = m_Clock->GetTimeStamp();
        double tillNow = (now - m_StartTime); // sec
        double leftTime = tillNow / prog * (100 - prog);

        if ((leftTime / 60.) >= 1.) // >= 1 min
        { 
          if ((leftTime / 3600.) >= 1.) // >= 1 hour
          { 
            if ((leftTime / 86400.) >= 1.) // >= 1 day
            { 
              double days = floor(leftTime / 86400.);
              double hours = floor((leftTime - days * 86400.) / 3600.);
              double min = floor((leftTime - days * 86400. -
                hours * 3600.) / 60.);
              double sec = leftTime - days * 86.400 - hours * 3600. -
                min * 60.;

              leftTimeStr.str(""); // clear
              leftTimeStr << days << "d " << hours << "h " <<
                min << "m " << itk::Math::Round<int, double>(sec) << "s";
            }
            else
            { 
              double hours = floor(leftTime / 3600.);
              double min = floor((leftTime - hours * 3600.) / 60.);
              double sec = leftTime - hours * 3600. - min * 60.;

              leftTimeStr.str(""); // clear
              leftTimeStr << hours << "h " << min << "m " <<
                itk::Math::Round<int, double>(sec) << "s";
            }
          }
          else
          { 
            double min = floor(leftTime / 60.);
            double sec = leftTime - min * 60.;

            leftTimeStr.str(""); // clear
            leftTimeStr << min << "m " << itk::Math::Round<int, double>(sec) << "s";
          }
        }
        else
        { 
          leftTimeStr.str(""); // clear
          leftTimeStr << leftTime << "s";
        }
      }
      std::cout << "\nPROGRESS: " << prog <<
        "% of iterations done.\t\t\tLEFT (for level " <<
        (m_Framework->GetRegistration()->GetCurrentLevel() + 1) <<
        "): " << leftTimeStr.str() << std::endl << std::endl;
    }
  }
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::WriteMovingImage(int currIt, bool costValueChanged)
{
  // projection-image-output each n-th iteration:
  if ((m_ImageOutput && (m_ImageModulo > 0) &&
       (currIt > 0) && m_Framework &&
       (m_ImageBaseFileName.length() > 0) &&
       ((currIt % m_ImageModulo) == 0)) ||
      (m_ImageOutput && m_ImageAutoOutput &&
       (currIt > 0) && m_Framework &&
       (m_ImageBaseFileName.length() > 0) &&
       costValueChanged))
  {
    if (m_Verbose)
      std::cout << std::endl << "Writing current moving image to disc (" <<
        "iteration=" << currIt << ")." << std::endl << std::endl;

    typedef typename FrameworkType::InternalImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    char imgfn[2048];
    std::string s = m_ImageBaseFileName;
    std::ostringstream ss("");

    ss << currIt; // insert iteration-number right before last "."
    std::string::size_type last = s.rfind('.');
    if (last != std::string::npos && last > 1)
      s.insert(last, ss.str());

    sprintf(imgfn, s.c_str(),
      (m_Framework->GetRegistration()->GetCurrentLevel() + 1));

    writer->SetFileName(imgfn);

    typename ImageType::ConstPointer actMovImg = m_Framework->
      GetCurrentProjectionImage();

    writer->SetInput(actMovImg);
    try
    { 
      writer->Update();
    }
    catch (ExceptionObject &e)
    { 
      std::cerr << "Error during writing current moving image (" <<
        imgfn << "): " << e << std::endl;
    }
  }
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::WritePDFImage(int currIt, bool costValueChanged)
{
  // metric-PDF-output each n-th iteration:
  if ((m_PDFOutput && (m_PDFModulo > 0) &&
       (currIt > 1) && m_Framework &&
       (m_PDFBaseFileName.length() > 0) &&
       ((currIt % m_PDFModulo) == 0)) ||
      (m_PDFOutput && m_PDFAutoOutput &&
       (currIt > 0) && m_Framework &&
       (m_PDFBaseFileName.length() > 0) &&
       costValueChanged))
  {
    if (m_Verbose)
      std::cout << std::endl << "Writing current PDF (metric) to disc (" <<
        "iteration=" << currIt << ")." << std::endl << std::endl;

    typedef typename FrameworkType::HistImageType ImageType;
    typedef ImageFileWriter<ImageType> WriterType;

    typename WriterType::Pointer writer = WriterType::New();

    char imgfn[2048];
    std::string s = m_PDFBaseFileName;
    std::ostringstream ss("");

    ss << currIt; // insert iteration-number right before last "."
    std::string::size_type last = s.rfind('.');
    if (last != std::string::npos && last > 1)
      s.insert(last, ss.str());

    sprintf(imgfn, s.c_str(),
      (m_Framework->GetRegistration()->GetCurrentLevel() + 1));

    writer->SetFileName(imgfn);

    typename ImageType::Pointer actPDFImg = m_Framework->
      GetCurrentMetricHistogram(false); // metric is up-to-date!

    if (actPDFImg) // could be NULL in case of non-histogram-based metrics
    { 
      actPDFImg->DisconnectPipeline();
      writer->SetInput(actPDFImg);
      try
      { 
        writer->Update();
      }
      catch (ExceptionObject &e)
      { 
        std::cerr << "Error during writing current histogram (" <<
          imgfn << "): " << e << std::endl;
      }
    }
  }
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::OpenCloseSplitLogFile(int currIt)
{
  std::ostringstream os; // string stream helper

  if (currIt == 1)
  { 
    if ((m_Framework->GetRegistration()->GetCurrentLevel() > 0) && m_Logging)
      m_LogFile.close();

    char fn[2048];
    os.str("");
    os << m_LogFileName;

    m_CurrLogSplitId++;
    if (m_LogSplitValue > 0)
      os << "." << m_CurrLogSplitId;

    sprintf(fn, os.str().c_str(),
      (m_Framework->GetRegistration()->GetCurrentLevel() + 1));
    m_LogFile.open(fn);
  }
  else if ((m_LogSplitValue > 0) && ((currIt % m_LogSplitValue) == 0))
  { 
    char fn[2048];

    m_CurrLogSplitId++;
    os.str("");
    os << m_LogFileName << "." << m_CurrLogSplitId;

    m_LogFile.close();

    // in case that log file name contains a "%d" for generic level string
    sprintf(fn, os.str().c_str(),
      (m_Framework->GetRegistration()->GetCurrentLevel() + 1));
    m_LogFile.open(fn);
  }
}

template <typename TFramework>
void 
OptimizerIterationCommand<TFramework>
::VerboseAndLogRegistrationParameters(int currIt,
  typename BaseMetricType::MeasureType thisIterationCostValue,
  typename BaseOptimizerType::ParametersType optimizerPos)
{
  std::ostringstream os; // string stream helper

  os.str(""); // clear

  // current iteration number:
  os << currIt << "\t";
  os << thisIterationCostValue;

  // current transform parameters:
  for (unsigned int i = 0; i < optimizerPos.GetSize(); ++i)
    os << "\t" << optimizerPos[i];

  if (m_LogCurrentParametersAlso) // apply additional logging parameters
  { 
    BaseTransformPointer transform = m_Framework->GetTransformation();
    typename BaseTransformType::ParametersType transPos =
      transform->GetParameters();
    BaseMetricPointer metric = m_Framework->GetMetric();
    typename BaseMetricType::MeasureType currCost =
      metric->GetValue(transPos); // this makes registration process slow ;(

    // append the additional parameters to log line:
    os << "\t" << currCost;
    for (unsigned int i = 0; i < transPos.GetSize(); ++i)
      os << "\t" << transPos[i];
  }

  if (m_Verbose)
    // logging-information also on stdout:
    std::cout << os.str() << std::endl;

  if (m_Logging && (m_LogFileName.length() > 0))
    m_LogFile << os.str() << std::endl;
}

template <typename TFramework>
bool 
OptimizerIterationCommand<TFramework>
::GetCurrentIterationParameters(BaseOptimizerConstPointer optimizer,
  int &currIt, int &maxIt,
  typename BaseMetricType::MeasureType &thisIterationCostValue,
  typename BaseOptimizerType::ParametersType &pos)
{
  bool validResults = false;

  if (!optimizer)
    return false;

  /* typedefs for concrete casting */
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

  if (optimizer->GetNameOfClass() == exhOptimizer->GetNameOfClass())
  { 
    exhOptimizer = SmartPointer<const ExhOptimizerType>(
      dynamic_cast<const ExhOptimizerType *>(optimizer.GetPointer()));
    currIt = exhOptimizer->GetCurrentIteration() + 1; // zero-based
    maxIt = exhOptimizer->GetMaximumNumberOfIterations();
    validResults = true;
    thisIterationCostValue = exhOptimizer->GetCurrentValue();
    pos = optimizer->GetCurrentPosition();
  }
  else if (optimizer->GetNameOfClass() == evolOptimizer->GetNameOfClass())
  { 
    evolOptimizer = SmartPointer<const EvolOptimizerType>(
      dynamic_cast<const EvolOptimizerType *>(optimizer.GetPointer()));
    currIt = evolOptimizer->GetCurrentIteration();
    maxIt = evolOptimizer->GetMaximumIteration();
    validResults = true;
    thisIterationCostValue = evolOptimizer->GetCurrentCost();
    pos = optimizer->GetCurrentPosition();
  }

  return validResults;
}

template <typename TFramework>
void
OptimizerIterationCommand<TFramework>
::StopOptimizer()
{ 
  if (!m_Framework)
    return;
  
  /* typedefs for concrete casting */
  typedef ExhaustiveOptimizerComplete ExhOptimizerType;
  typedef ExhOptimizerType::Pointer ExhOptimizerPointer;
  typedef OnePlusOneEvolutionaryOptimizer EvolOptimizerType;
  typedef EvolOptimizerType::Pointer EvolOptimizerPointer;
  
  ExhOptimizerPointer exh = ExhOptimizerType::New();
  EvolOptimizerPointer evol = EvolOptimizerType::New();
  
  BaseOptimizerPointer opt = m_Framework->GetOptimizer();
  
  if (!opt)
    return;
  
  if (opt->GetNameOfClass() == exh->GetNameOfClass())
  {
    exh = SmartPointer<ExhOptimizerType>(
      dynamic_cast<ExhOptimizerType *>(opt.GetPointer()));
    
    exh->StopWalking();
  }
  else if (opt->GetNameOfClass() == evol->GetNameOfClass())
  {
    evol = SmartPointer<EvolOptimizerType>(
      dynamic_cast<EvolOptimizerType *>(opt.GetPointer()));
      
    evol->StopOptimization();
  }
}


}


#endif /* ITKOPTIMIZERITERATIONCOMMAND_TXX_ */
