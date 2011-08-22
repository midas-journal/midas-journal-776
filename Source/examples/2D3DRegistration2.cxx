//
#include <iostream>
#include <fstream>
#include <string>

#include <itksys/SystemTools.hxx>
#include <itkImageFileReader.h>

#include <itkRayCastPerspectiveProjectionImageFilter.h>
#include <itkSplatPerspectiveProjectionImageFilter.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>
#include <itkOnePlusOneEvolutionaryOptimizer.h>
#include <itkExhaustiveOptimizerComplete.hxx>
#include <itkMattesMutualInformationImageToImageMetricComplete.h>
#include <itkRegularStepGradientDescentOptimizer.h>

#include "oraMultiResolutionImage2D3DRegistrationMethodWithSRC.h"
#include "oraStochasticRankCorrelationImageToImageMetric.h"

using namespace itk;


/**
  * Display usage information (command line parameter options) and some
  * additional remarks.
  */
void Usage()
{
  std::cout << "*********************************************" << std::endl;
  std::cout << "***     2D/3D-REGISTRATION (WITH SRC)     ***" << std::endl;
  std::cout << "*********************************************" << std::endl;

  std::cout << std::endl;
  std::cout << "Author: Philipp Steininger" << std::endl;
  std::cout << "Affiliation: radART Institute, PMU, Salzburg, Austria" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Application usage: 2D3DRegistration2 <config-file> " <<
    "<fixed-image> <input-volume>" << std::endl;
  std::cout << std::endl;
  std::cout << "<config-file> is the configuration file containing " <<
    "different keys and values controlling the program (have a look at " <<
    "\"RegistrationExample.cfg\")." << std::endl;
  std::cout << "<fixed-image> is the reference 2D image (e.g. an X-ray); " <<
    "metaheader-format (*.mhd) preferred" << std::endl;
  std::cout << "<input-volume> is the (static) 3D image (e.g. a CT image); " <<
    "metaheader-format (*.mhd) preferred" << std::endl;
  std::cout << std::endl;
}

/**
 * Load the specified image file with registration framework's internal image
 * type (float, 3D).
 * @param fileName file name of image to be loaded
 * @return pointer to loaded image if successful, NULL otherwise
 */
MultiResolutionImage2D3DRegistrationMethod<>::InternalImageType::Pointer
LoadImage(std::string fileName)
{
  // Image type
  typedef MultiResolutionImage2D3DRegistrationMethod<>::InternalImageType
    ImageType;
  ImageType::Pointer image;

  try // try to load the image using the ITK-reader-factory
  {
    ImageFileReader<ImageType>::Pointer reader =
      ImageFileReader<ImageType>::New();

    reader->SetFileName(fileName.c_str());
    reader->Update();

    image = reader->GetOutput();
    image->DisconnectPipeline(); // disconnect to avoid inconsistencies later
  }
  catch (ExceptionObject &e) // error
  {
    image = NULL;
  }

  return image;
}

/**
 * Look if the specified path exists. If it does not exist the procedure creates
 * it.
 * @param pathName path specification - should really be a path without file
 * specification!
 * @return whether path exists now (true if it could be created or had already
 * existed)
 */
bool EnsurePathExists(std::string pathName)
{
  if (!itksys::SystemTools::FileExists(pathName.c_str(), false)) // does not ex.
  {
    itksys::SystemTools::MakeDirectory(pathName.c_str()); // create tree

    return itksys::SystemTools::FileExists(pathName.c_str(), false); // check
  }
  else // had already existed
    return true;
}

/**
 * Read a specified key from the specified configuration file.
 * Format: key=value (case sensitive and without leading spaces!)
 * @param key specifies the key (case sensitive)
 * @param configFile path and file name of the configuration file
 * @return value of the key if found or empty string if key could not be found
 */
std::string GetKeyFromConfigFile(const std::string key,
  const std::string configFile)
{
  std::ifstream config;
  char linebuff[2048];

  config.open(configFile.c_str()); // open file

  if (config) // success
  {
    std::string resultString = "";

    while (config.getline(linebuff, 2048) &&
           resultString.length() == 0) // read line-by-line
    {
      std::string line(linebuff); // convert into string

      if (line.find(key + "=") == 0) // search for separating "="
        resultString = line.substr(key.length() + 1, line.length()); // cut
    }

    config.close(); // ready

    return resultString;
  }
  else // file could not be opened - error
    return "";
}


/**
 * <b>2D3DRegistration.</b> <br>
 * A 2D/3D-registration application which aims at aligning a 3D static volume
 * with a 2D image outgoing from a manually configured initial position. <br>
 * The application is configured and parametrized by a configuration file (
 * have a look at the exemplary file "RegistrationExample.cfg" in the source
 * folder).
 * @param argc number of arguments (argv); at least 4 arguments are expected
 * @param argv array of arguments; the first argument is automatically applied
 * (name of the application); the second argument is expected to be the file
 * name of the configuration file; the third argument is expected to be the
 * file name of the fixed (2D) image; the fourth argument is expected to be the
 * file name of the moving (3D) image
 */
int main(int argc, char *argv[])
{
  if (argc < 4) // check arguments - show usage-information on failure
  {
    Usage();
    return EXIT_FAILURE;
  }

  // store the argument values:
  std::string configFile = argv[1];
  std::string fixedImageFile = argv[2];
  std::string movingImageFile = argv[3];


  // create the registration framework (with default template parameters):
  typedef ora::MultiResolutionImage2D3DRegistrationMethodWithSRC<> RegistrationType;
  typedef RegistrationType::Pointer RegistrationPointer;

  RegistrationPointer registration = RegistrationType::New();


  // load and set the images:
  typedef ora::MultiResolutionImage2D3DRegistrationMethodWithSRC<>::InternalImageType
    ImageType;
  typedef ImageType::Pointer ImagePointer;
  typedef ora::MultiResolutionImage2D3DRegistrationMethodWithSRC<>::RankImageType
    RankImageType;

  ImagePointer fixedImage = LoadImage(fixedImageFile); // 2D
  ImagePointer movingImage = LoadImage(movingImageFile); // 3D

  registration->SetFixedImage(fixedImage);
  registration->SetMovingImage(movingImage);


  // set the image regions according to the configuration:
  typedef RegistrationType::InternalRegionType RegionType;
  typedef RegistrationType::InternalSpacingType SpacingType;

  SpacingType xspacing; // spacing of the DRR
  RegionType xregion; // region of the DRR projection

  xspacing[0] = atof(GetKeyFromConfigFile("xrayspacing[0]",
    configFile).c_str());
  xspacing[1] = atof(GetKeyFromConfigFile("xrayspacing[1]",
    configFile).c_str());
  for (int d = 0; d < 2; d++) // apply default spacing on invalid values
    if (xspacing[d] < 1e-6)
      xspacing[d] = 1.;
  xregion.SetIndex(0, itk::Math::Round<int, double>(atof(GetKeyFromConfigFile(
      "xrayplaneoffset_mm[0]", configFile).c_str()) / xspacing[0]));
  xregion.SetIndex(1, itk::Math::Round<int, double>(atof(GetKeyFromConfigFile(
      "xrayplaneoffset_mm[1]", configFile).c_str()) / xspacing[1]));
  xregion.SetSize(0, itk::Math::Round<int, double>(atof(GetKeyFromConfigFile("xrayplanesize_mm[0]",
    configFile).c_str()) / xspacing[0]));
  xregion.SetSize(1, itk::Math::Round<int, double>(atof(GetKeyFromConfigFile("xrayplanesize_mm[1]",
    configFile).c_str()) / xspacing[1]));

  registration->SetXrayImageRegion(xregion);
  registration->SetXrayImageSpacing(xspacing);
  // take the X-ray image region over as fixed image region of registration
  registration->SetFixedImageRegion(xregion);


  // instantiate the basic registration components:
  // component-selections from config-file:
  std::string projectortype = GetKeyFromConfigFile("projectortype", configFile);
  std::string metrictype = GetKeyFromConfigFile("metrictype[0]", configFile);
  std::string optimizertype = GetKeyFromConfigFile("optimizertype[0]",
    configFile);
  std::string transformtype = GetKeyFromConfigFile("transformtype", configFile);

  // - projector component:
  RegistrationType::BaseProjectionFilterPointer projector = NULL;

  if (projectortype == "RAYC") // ray-casting
    projector = static_cast<RegistrationType::BaseProjectionFilterType *>(
      RayCastPerspectiveProjectionImageFilter<ImageType, ImageType, double,
      Collector::SumCollector<ImageType::PixelType> >::New().GetPointer());
  else if (projectortype == "SPLAT") // wobbled splatting
    projector = static_cast<RegistrationType::BaseProjectionFilterType *>(
      SplatPerspectiveProjectionImageFilter<ImageType, ImageType, double,
      Collector::SumCollector<ImageType::PixelType> >::New().GetPointer());

  registration->SetProjectionFilter(projector);

  // - metric component:
  RegistrationType::BaseMetricPointer metric = NULL;

  if (metrictype == "NMI") // normalized mutual information
    metric = static_cast<RegistrationType::BaseMetricType *>(
      NormalizedMutualInformationHistogramImageToImageMetric<ImageType,
      ImageType>::New().GetPointer());
  else if (metrictype == "MMI") // Mattes mutual information
    metric = static_cast<RegistrationType::BaseMetricType *>(
      MattesMutualInformationImageToImageMetricComplete<ImageType,
      ImageType>::New().GetPointer());
  else if (metrictype == "SRC") // stochastic rank correlation
  {
    ora::StochasticRankCorrelationImageToImageMetric<ImageType,
          ImageType, RankImageType>::Pointer src =
              ora::StochasticRankCorrelationImageToImageMetric<ImageType,
                    ImageType, RankImageType>::New();
    if (atoi(GetKeyFromConfigFile("pdfoutput", configFile).c_str()))
      src->SetExtractSampleDistribution(true); // debugging
    metric = static_cast<RegistrationType::BaseMetricType *>(
      src.GetPointer());
  }

  registration->SetMetric(metric);

  // - optimizer component:
  RegistrationType::BaseOptimizerPointer optimizer = NULL;

  if (optimizertype == "EXH") // exhaustive search space exploration
    optimizer = static_cast<RegistrationType::BaseOptimizerType *>(
      ExhaustiveOptimizerComplete::New().GetPointer());
  else if (optimizertype == "EVOL") // evolutionary 1+1 optimization
    optimizer = static_cast<RegistrationType::BaseOptimizerType *>(
      OnePlusOneEvolutionaryOptimizer::New().GetPointer());
  else if (optimizertype == "RSGD") // regular step gradient descent opt.
    optimizer = static_cast<RegistrationType::BaseOptimizerType *>(
      RegularStepGradientDescentOptimizer::New().GetPointer());

  registration->SetOptimizer(optimizer);

  // - transformation component:
  RegistrationType::BaseTransformPointer transformation = NULL;

  if (transformtype == "E3D") // Euler 3D transformation (rigid)
    transformation = static_cast<RegistrationType::BaseTransformType *>(
      Euler3DTransform<RegistrationType::ScalarType>::New().GetPointer());
  else if (transformtype == "AFF") // affine transformation
    transformation = static_cast<RegistrationType::BaseTransformType *>(
      AffineTransform<RegistrationType::ScalarType>::New().GetPointer());

  registration->SetTransformation(transformation);


  registration->Initialize(); // tie the registration components together


  bool success = true; // success-state of configure-methods

  // configure the framework now by applying generic strings:

  // - set up the images:
  std::string imagesconfig = GetKeyFromConfigFile("imagesconfig", configFile);
  // configure string FORMAT: "volumeorigin <double> <double> <double>
  // fixedrescaling <double> <double>"
  success &= registration->ConfigureImages(imagesconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize images (" <<
      imagesconfig << ")!" << std::endl;

  // - set up the image pyramids:
  std::string pyramidsconfig = GetKeyFromConfigFile("pyramidsconfig",
    configFile);
  // configure string format: "levels <int>"
  success &= registration->ConfigureImagePyramids(pyramidsconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize pyramids (" <<
      pyramidsconfig << ")!" << std::endl;

  // - set up the transformation:
  std::string transformationconfig = GetKeyFromConfigFile(
    "transformationconfig", configFile);
  // configure string format:
  // (E3D) "cor <double> <double> <double> rot <double> <double> <double>
  // transl <double> <double> <double>"
  // (AFF) "cor <double> <double> <double> rot <double> <double> <double>
  // transl <double> <double> <double> scale <double> <double> <double>"
  success &= registration->ConfigureTransformation(transformationconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize transformation (" <<
      transformationconfig << ")!" << std::endl;

  // - set up the optimizer:
  std::string optimizerconfig = GetKeyFromConfigFile("optimizerconfig[0]",
    configFile);
  // configure string format:
  // (EXHAUSTIVE)
  // "min <bool> scales <double> <double> ... steplen <double>
  // stepcount <int> <int> ..."
  // (EVOLUTIONARY)
  // "min <bool> scales <double> <double> ... iterations <int> nvgseed <int>
  // radius <double> gfact <double> sfact <double> epsilon <double>"
  // (REGSTEPGRADDESC)
  // "min <bool> scales <double> ... <double> steps <double> <double>
  // iterations <int> gradtol <double> relax <double>"
  success &= registration->ConfigureOptimizer(optimizerconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize optimizer (" <<
      optimizerconfig << ")!" << std::endl;

  // - set up the metric:
  std::string metricconfig = GetKeyFromConfigFile("metricconfig[0]",
    configFile);
  // configure string format:
  // (NMI)
  // "histbins <int> <int> lowerbounds <double> <double>
  // upperbounds <double> <double>"
  // (MMI)
  // "histbins <int> lowerbounds <double> <double>
  // upperbounds <double> <double> nosamples <double>"
  // (RSGD)
  // (SRC)
  // "rseeds <int> <int> dscales <double> ... <double>
  // fhist <double> <double> <int> <bool> mhist <double> <double> <int> <bool>
  // coverage <double>"
  success &= registration->ConfigureMetric(metricconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize metric (" <<
      metricconfig << ")!" << std::endl;

  // - set up the interpolator (projection-component):
  std::string interpolatorconfig = GetKeyFromConfigFile("interpolatorconfig[0]",
    configFile);
  // configure string format:
  // (SPLATTING)
  // "jobs <int> threshold <double> transferfunc <int> <double> <double> ...
  // transferinterp <int> focalpoint <double> <double> <double>
  // rescale <double> <double> wobbmode <int> wobbsigma <double>
  // wobbcount <int> splatpixassign <int> wobbdeterm <bool>
  // rescale <double> <double>"
  // (RAYCASTING)
  // "jobs <int> threshold <double> transferfunc <int> <double> <double> ...
  // transferinterp <int> focalpoint <double> <double> <double>
  // rescale <double> <double> raystepsize <double> raystepinterp <int>"
  success &= registration->ConfigureInterpolator(interpolatorconfig);
  if (!success)
    std::cerr << "CONFIGURE-ERROR: Could not initialize interpolator (" <<
      interpolatorconfig << ")!" << std::endl;

  // compute the necessary translations resulting from ref.-point-projection:
  // configure string format:

  if (atoi(GetKeyFromConfigFile("usereferencepoint", configFile).c_str()))
  {
    std::string referencepointconfig = GetKeyFromConfigFile(
      "referencepointconfig", configFile);
    // "inplane <double> <double> refpoint <double> <double> <double>
    // applytransl <bool>"
    success &= registration->ComputeVolumeTranslationByReferencePoint(
      referencepointconfig);
    if (!success)
      std::cerr << "CONFIGURE-ERROR: Could not compute volume translation (" <<
        referencepointconfig << ")!" << std::endl;
  }

  // now check if all configuration steps have been successful, hook in the
  // observers and start the registration:
  if (success)
  {
    // hook in observers (registration framework, multiresolution, optimizer):
    typedef RegistrationType::GenericFrameworkObserverType
      FrameworkObserverType;
    typedef RegistrationType::MultiResolutionObserverType RegObserverType;
    typedef RegistrationType::OptimizationObserverType OptObserverType;

    // basic output path for documentation of the registration:
    std::string outputPath = GetKeyFromConfigFile("outputbasepath", configFile);
    EnsurePathExists(outputPath); // create path if necessary

    // - registration framework observer:
    FrameworkObserverType::Pointer fwObserver = FrameworkObserverType::New();

    fwObserver->SetVerbose(atoi(GetKeyFromConfigFile("verbose",
      configFile).c_str())); // verbose related messages
    fwObserver->SetLogging(atoi(GetKeyFromConfigFile("logging",
      configFile).c_str())); // log to a log-file
    fwObserver->SetLogFileName(outputPath +
      GetKeyFromConfigFile("resultlogfile", configFile)); // final (last level)
    fwObserver->SetImageOutput(atoi(GetKeyFromConfigFile(
      "imageoutput", configFile).c_str()));
    std::vector<std::string> fnames;
    fnames.push_back(outputPath + GetKeyFromConfigFile("resultmovingimagefile",
      configFile)); // moving image file names
    fnames.push_back(outputPath + GetKeyFromConfigFile("resultmovingpdffile",
      configFile)); // pdf image file names
    fwObserver->SetImageFileNames(fnames); // take over

    registration->AddGenericFrameworkObserver(fwObserver); // add the observer

    // - multiresolution observer:
    RegObserverType::Pointer regObserver = RegObserverType::New();

    regObserver->SetVerbose(atoi(GetKeyFromConfigFile("verbose",
      configFile).c_str())); // verbose related messages
    regObserver->SetLogging(atoi(GetKeyFromConfigFile("logging",
      configFile).c_str())); // log to a log-file
    regObserver->SetLogFileName(outputPath + GetKeyFromConfigFile(
      "resultlevellogfile", configFile)); // level results (except last level)

    fnames.clear(); // image output file names:
    fnames.push_back(outputPath + GetKeyFromConfigFile("fixedimagefile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("inputvolumefile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("intialmovingimagefile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("finalmovingimagefile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("finalmovingpdffile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("croppedfixedimagefile",
      configFile));
    fnames.push_back(outputPath + GetKeyFromConfigFile("fixedimagemaskfile",
      configFile));
    regObserver->SetImageOutput(atoi(GetKeyFromConfigFile(
      "imageoutput", configFile).c_str()));
    regObserver->SetImageFileNames(fnames); // take over
    regObserver->SetFramework(registration); // set reference to framework

    std::vector<std::string> mcfgs;
    std::vector<std::string> msels;
    int levels = registration->GetRegistration()->GetNumberOfLevels();
    char buff[100];
    for (int i = 1; i < levels; i++)
    {
      sprintf(buff, "metrictype[%d]", i);
      std::string cfg = GetKeyFromConfigFile(buff, configFile);
      if (cfg.length() > 0) // require consecutive entries!
        msels.push_back(cfg);
      else
        break;
    }
    regObserver->SetMetricSels(msels);
    for (int i = 1; i < levels; i++)
    {
      sprintf(buff, "metricconfig[%d]", i);
      std::string cfg = GetKeyFromConfigFile(buff, configFile);
      if (cfg.length() > 0) // require consecutive entries!
        mcfgs.push_back(cfg);
      else
        break;
    }
    regObserver->SetMetricConfigs(mcfgs);
    std::vector<std::string> ocfgs;
    std::vector<std::string> osels;
    for (int i = 1; i < levels; i++)
    {
      sprintf(buff, "optimizertype[%d]", i);
      std::string cfg = GetKeyFromConfigFile(buff, configFile);
      if (cfg.length() > 0) // require consecutive entries!
        osels.push_back(cfg);
      else
        break;
    }
    regObserver->SetOptimizerSels(osels);
    for (int i = 1; i < levels; i++)
    {
      sprintf(buff, "optimizerconfig[%d]", i);
      std::string cfg = GetKeyFromConfigFile(buff, configFile);
      if (cfg.length() > 0) // require consecutive entries!
        ocfgs.push_back(cfg);
      else
        break;
    }
    regObserver->SetOptimizerConfigs(ocfgs);
    std::vector<std::string> icfgs;
    for (int i = 1; i < levels; i++)
    {
      sprintf(buff, "interpolatorconfig[%d]", i);
      std::string cfg = GetKeyFromConfigFile(buff, configFile);
      if (cfg.length() > 0) // require consecutive entries!
        icfgs.push_back(cfg);
      else
        break;
    }
    regObserver->SetInterpolatorConfigs(icfgs);

    registration->AddMultiResolutionObserver(regObserver); // add the observer

    // - optimization observer:
    OptObserverType::Pointer optObserver = OptObserverType::New();

    optObserver->SetVerbose(atoi(GetKeyFromConfigFile("verbose",
      configFile).c_str())); // verbose related messages
    optObserver->SetLogging(atoi(GetKeyFromConfigFile("logging",
      configFile).c_str())); // log to a log-file
    optObserver->SetLogCurrentParametersAlso(atoi(GetKeyFromConfigFile(
      "logcurrentparametersalso", configFile).c_str())); // actual trans. par.
    optObserver->SetLogFileName(outputPath + GetKeyFromConfigFile(
      "iterationparameterslogfile", configFile).c_str()); // iteration log file
    optObserver->SetLogSplitValue(atoi(GetKeyFromConfigFile(
      "logsplitvalue", configFile).c_str())); // split log file or not
    optObserver->SetImageOutput(atoi(GetKeyFromConfigFile(
      "imageoutput", configFile).c_str())); // image output
    optObserver->SetImageModulo(atoi(GetKeyFromConfigFile(
      "imagemodulo", configFile).c_str())); // each n-th image is written
    optObserver->SetImageAutoOutput(atoi(GetKeyFromConfigFile(
      "imageautooutput", configFile).c_str())); // automatic detection
    optObserver->SetImageBaseFileName(outputPath +
      GetKeyFromConfigFile("actualmovingimagefile", configFile)); // mov. img.
    optObserver->SetPDFOutput(atoi(GetKeyFromConfigFile(
      "pdfoutput", configFile).c_str())); // PDF output
    optObserver->SetPDFModulo(atoi(GetKeyFromConfigFile(
      "pdfmodulo", configFile).c_str())); // each n-th PDF is written
    optObserver->SetPDFAutoOutput(atoi(GetKeyFromConfigFile(
      "pdfautooutput", configFile).c_str())); // automatic detection
    optObserver->SetPDFBaseFileName(outputPath +
      GetKeyFromConfigFile("actualpdfimagefile", configFile)); // PDF file
    optObserver->SetFramework(registration); // set reference to framework

    registration->AddOptimizationObserver(optObserver); // add the observer


    // execute the registration (as filter)
    try
    {
      std::cout << std::endl << "starting registration ..." << std::endl;
      registration->Update();
      std::cout << std::endl << "... finished registration" << std::endl;
    }
    catch(ExceptionObject &err)
    {
      std::cerr << "Error during Registration: " << std::endl;
      std::cerr << err << std::endl;

      return EXIT_FAILURE;
    }

    std::cout << "Successfully finished the registration." << std::endl;

    return EXIT_SUCCESS;
  }
  else
  {
    std::cerr << "There have been too many errors, cannot start " <<
        "registration." << std::endl;

    return EXIT_FAILURE;
  }
}
