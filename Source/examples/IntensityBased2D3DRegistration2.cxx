
// NOTE:
//
// This simple intensity-based 2D/3D-registration application is a modified
// version of "IntensityBased2D3DRegistration" from InsightApplications.
//
// It was originally written by
//  john.hipwell@kcl.ac.uk
//   and
//  thomas@hartkens.de.
//
// This application utilizes the stochastic rank correlation (SRC) metric and
// regular step gradient descent optimization. For ray-casting,
// itk::RayCastInterpolateImageFunction is used.
//
// Have a look at the command line help; there are new command line arguments!
//
// Philipp Steininger
// Institute for Research and Development on Advanced Radiation Technologies (radART)
// Paracelsus Medical University (PMU)
// Salzburg, Austria
//
// Markus Neuner
// Institute for Research and Development on Advanced Radiation Technologies (radART)
// Paracelsus Medical University (PMU)
// Salzburg, Austria


#include <iostream>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkEuler3DTransform.h>
#include "oraStochasticRankCorrelationImageToImageMetric.h"
#include <itkRayCastInterpolateImageFunction.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkResampleImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkSquaredDifferenceImageFilter.h>
#include <itkCommand.h>
#include <itkImageMaskSpatialObject.h>
#include <itkMaskImageFilter.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


// Registration Monitor:
class CommandIterationUpdate : public itk::Command 
{
public:
  typedef  CommandIterationUpdate   Self;
  typedef  itk::Command             Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() {};

public:
  typedef itk::RegularStepGradientDescentOptimizer     OptimizerType;
  typedef const OptimizerType                         *OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object *object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = 
                      dynamic_cast< OptimizerPointer >( object );
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
        return;
      }
    std::cout << "Iteration: " << optimizer->GetCurrentIteration() << std::endl;
    std::cout << " Similarity: " << optimizer->GetValue() << std::endl;
    std::cout << " Position: " << optimizer->GetCurrentPosition() << std::endl;
    std::cout << " Current Step Size: " << optimizer->GetCurrentStepLength() << std::endl;
  }
};

// Multi-resolution monitor:
double StepLengthAdjustmentFactor = 0.5; // multiplication factor
template <typename TRegistration>
class RegistrationInterfaceCommand : public itk::Command 
{
public:
  typedef  RegistrationInterfaceCommand   Self;
  typedef  itk::Command                   Superclass;
  typedef  itk::SmartPointer<Self>        Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand() {};

public:
  typedef   TRegistration                              RegistrationType;
  typedef   RegistrationType *                         RegistrationPointer;
  typedef   itk::RegularStepGradientDescentOptimizer   OptimizerType;
  typedef   OptimizerType *                            OptimizerPointer;

  void Execute(itk::Object * object, const itk::EventObject & event)
  {
    if( typeid( event ) != typeid( itk::IterationEvent ) )
      {
      return;
      }

    RegistrationPointer registration =
                            dynamic_cast<RegistrationPointer>( object );

    // Step length adjustment:
    OptimizerPointer optimizer = dynamic_cast< OptimizerPointer >( 
                       registration->GetOptimizer() );

    if ( registration->GetCurrentLevel() != 0 )
      {
      optimizer->SetMaximumStepLength(optimizer->GetMaximumStepLength() * StepLengthAdjustmentFactor);
      optimizer->SetMinimumStepLength(optimizer->GetMinimumStepLength() * StepLengthAdjustmentFactor);
      }

    optimizer->Print( std::cout );
  }

  void Execute(const itk::Object * , const itk::EventObject & )
    { return; }
};

// Program usage information:
void usage()
{
  std::cerr << "\n";
  std::cerr << "Usage: itkMultiResIntensityBasedRegn2D3D <options> Image2D Volume3D\n";
  std::cerr << "       Registers a 3D volume to a 2D image. \n\n";
  std::cerr << "   where <options> is one or more of the following:\n\n";
  std::cerr << "       <-h>                     Display (this) usage information\n";
  std::cerr << "       <-v>                     Verbose output [default: no]\n";
  std::cerr << "       <-dbg>                   Debugging output [default: no]\n";
  std::cerr << "       <-n int>                 The number of scales to apply [default: 2]\n";
  std::cerr << "       <-maxScale int>          The scale factor corresponding to max resolution [default: 1]\n";
  std::cerr << "       <-step float float>      Maximum and minimum step sizes [default: 4 and 0.01]\n";
  std::cerr << "       <-fl float>              Focal length or source to image distance [default: 400mm]\n";
  std::cerr << "       <-t float float float>   Translation parameter of the camera \n";
  std::cerr << "       <-rx float>              Rotation around x,y,z axis in degrees \n";
  std::cerr << "       <-ry float>\n";
  std::cerr << "       <-rz float>\n";
  std::cerr << "       <-normal float float>    The 2D projection normal position [default: 0x0mm]\n";
  std::cerr << "       <-cor float float float> The centre of rotation relative to centre of volume\n";
  std::cerr << "       <-threshold float>       Threshold [default: 0]\n";
  std::cerr << "       <-o file>                Output image filename\n";
  std::cerr << "       <-r radius>              Fixed image mask with specified radius in PIXELS [default: 0, off<=0]\n";
  std::cerr << "       <-om file>               Output image filename masked [only performed if radius > 0]\n\n";
  std::cerr << "       <-initialpars float float float float float>  Initial transformation parameters directly [default: not specified]\n";
  std::cerr << "       <-maxit int> Maximum number of optimization iterations [default: 300]\n";
  std::cerr << "       <-rf float> Optimizer relaxation factor (step length adaption) between 0.0 and 1.0 [default: 0.5]\n";
  std::cerr << "       <-oscale float float float float float float> The optimizer scales for optimization [default: 180./M_PI 180./M_PI 180./M_PI 1 1 1]\n";
  std::cerr << "       <-samplecoverage float> The sample coverage for stochastic rank correlation [default: 10]\n";
  std::cerr << "       <-dscale float float float float float float> The derivative scales for gradient computation [default: 0.01745 0.01745 0.01745 1 1 1]\n";
  std::cerr << "       <-fixedhist float float int bool> The fixed histogram configuration: min-intensity, max-intensity, number of bins, clip at ends flag [default: 0 255 256 1]\n";
  std::cerr << "       <-movinghist float float int bool> The moving histogram configuration: min-intensity, max-intensity, number of bins, clip at ends flag [default: 0 255 256 1]\n";
  std::cerr << "       <-rseeds int int int> Seeds for deterministic stochastic sampling (mask) [default (non-deterministic): 0 0 0]\n";
  std::cerr << "       <-horn> use Horn-correction for tied ranks [default: do not use Horn-correction]\n";
  std::cerr << "       <-zeroRanksContribute> moving zero-ranks contribute to coefficient [default: no, they do not contribute]\n";
  std::cerr << "       <-io string> Initial projection image output (file name) [default: NULL]\n\n";
  std::cerr << "                                by  john.hipwell@kcl.ac.uk\n";
  std::cerr << "                                and thomas@hartkens.de\n";
  std::cerr << "                                (Imaging Sciences KCL London)\n\n";
  std::cerr << "                                and philipp steininger\n";
  std::cerr << "                                (radART, PMU, Salzburg, Austria)\n";
  std::cerr << "                                and markus neuner\n";
  std::cerr << "                                (radART, PMU, Salzburg, Austria)\n\n";
  exit(1);
}

// Program:
int main( int argc, char *argv[] )
{
  char *fileImage2D = NULL;
  char *fileVolume3D = NULL;
  char *fileOutput = NULL;
  char *fileOutputMasked = NULL;
  char *initialOutput = NULL;
  bool ok;
  bool verbose = false;
  bool debug = false;
  unsigned int nScales = 2;
  int maxScale = 1;
  double rx = 0.;
  double ry = 0.;
  double rz = 0.;
  double tx = 0.;
  double ty = 0.;
  double tz = 0.;
  double cx = 0.;
  double cy = 0.;
  double cz = 0.;
  double focalLength = 400.;
  double maxStepSize = 4.;
  double minStepSize = 1.;
  double o2Dx = 0;
  double o2Dy = 0;
  double threshold=0;
  double dscale[] = {0.01745329251994329577, 0.01745329251994329577, 0.01745329251994329577, 1, 1, 1};
  double fhistmin = 0;
  double fhistmax = 255;
  int fhistbins = 256;
  bool fhistclip = true;
  double mhistmin = 0;
  double mhistmax = 255;
  int mhistbins = 256;
  bool mhistclip = true;
  double samplecoverage = 10;
  double oscale[] = {M_PI / 180, M_PI / 180, M_PI / 180, 1, 1, 1};
  int maxit = 300;
  unsigned int rseeds[] = {0, 0, 0};
  double rf = 0.5;
  double *initpars = NULL;
  bool horn = false;
  bool zeroRanksContribute = false;
  double fmRadius = 0;

  // Parse command line parameters
  if (argc <= 1)
    usage();

  while (argc > 1) 
    {
    ok = false;

    if ((ok == false) && (strcmp(argv[1], "-h") == 0))
      {
      argc--; argv++;
      ok = true;
      usage();      
      }

    if ((ok == false) && (strcmp(argv[1], "-v") == 0))
      {
      argc--; argv++;
      ok = true;
      verbose = true;
      }

    if ((ok == false) && (strcmp(argv[1], "-dbg") == 0))
      {
      argc--; argv++;
      ok = true;
      debug = true;
      }

    if ((ok == false) && (strcmp(argv[1], "-n") == 0))
      {
      argc--; argv++;
      ok = true;
      nScales=atoi(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-maxScale") == 0))
      {
      argc--; argv++;
      ok = true;
      maxScale=atoi(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-fl") == 0))
      {
      argc--; argv++;
      ok = true;
      focalLength=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-step") == 0))
      {
      argc--; argv++;
      ok = true;
      maxStepSize=atof(argv[1]);
      argc--; argv++;
      minStepSize=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-t") == 0))
      {
      argc--; argv++;
      ok = true;
      tx=atof(argv[1]);
      argc--; argv++;
      ty=atof(argv[1]);
      argc--; argv++;
      tz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-rx") == 0))
      {
      argc--; argv++;
      ok = true;
      rx=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-ry") == 0))
      {
      argc--; argv++;
      ok = true;
      ry=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-rz") == 0))
      {
      argc--; argv++;
      ok = true;
      rz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-normal") == 0))
      {
      argc--; argv++;
      ok = true;
      o2Dx=atof(argv[1]);
      argc--; argv++;
      o2Dy=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-cor") == 0))
      {
      argc--; argv++;
      ok = true;
      cx=atof(argv[1]);
      argc--; argv++;
      cy=atof(argv[1]);
      argc--; argv++;
      cz=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-threshold") == 0))
      {
      argc--; argv++;
      ok = true;
      threshold=atof(argv[1]);
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-o") == 0))
      {
      argc--; argv++;
      ok = true;
      fileOutput = argv[1];
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-om") == 0))
      {
      argc--; argv++;
      ok = true;
      fileOutputMasked = argv[1];
      argc--; argv++;
      }

    if ((ok == false) && (strcmp(argv[1], "-io") == 0))
      {
      argc--; argv++;
      ok = true;
      initialOutput = argv[1];
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-dscale") == 0))
      {
      argc--; argv++;
      ok = true;
      for (int xx = 0; xx < 6; xx++)
        {
        dscale[xx]=atof(argv[1]);
        argc--; argv++;
        }
      }
    if ((ok == false) && (strcmp(argv[1], "-fixedhist") == 0))
      {
      argc--; argv++;
      ok = true;
      fhistmin=atof(argv[1]);
      argc--; argv++;
      fhistmax=atof(argv[1]);
      argc--; argv++;
      fhistbins=atoi(argv[1]);
      argc--; argv++;
      fhistclip=(bool)atoi(argv[1]);
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-movinghist") == 0))
      {
      argc--; argv++;
      ok = true;
      mhistmin=atof(argv[1]);
      argc--; argv++;
      mhistmax=atof(argv[1]);
      argc--; argv++;
      mhistbins=atoi(argv[1]);
      argc--; argv++;
      mhistclip=(bool)atoi(argv[1]);
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-samplecoverage") == 0))
      {
      argc--; argv++;
      ok = true;
      samplecoverage=atof(argv[1]);
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-oscale") == 0))
      {
      argc--; argv++;
      ok = true;
      for (int xx = 0; xx < 6; xx++)
        {
        oscale[xx]=atof(argv[1]);
        argc--; argv++;
        }
      }
    if ((ok == false) && (strcmp(argv[1], "-initialpars") == 0))
      {
      argc--; argv++;
      ok = true;
      initpars = new double[6];
      for (int xx = 0; xx < 6; xx++)
        {
        initpars[xx]=atof(argv[1]);
        argc--; argv++;
        }
      }
    if ((ok == false) && (strcmp(argv[1], "-maxit") == 0))
      {
      argc--; argv++;
      ok = true;
      maxit=atoi(argv[1]);
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-rf") == 0))
      {
      argc--; argv++;
      ok = true;
      rf=atof(argv[1]);
      argc--; argv++;
      }
    if ((ok == false) && (strcmp(argv[1], "-rseeds") == 0))
      {
      argc--; argv++;
      ok = true;
      for (int xx = 0; xx < 3; xx++)
        {
        rseeds[xx]=atoi(argv[1]);
        argc--; argv++;
        }
      }
    if ((ok == false) && (strcmp(argv[1], "-horn") == 0))
      {
      argc--; argv++;
      ok = true;
      horn = true;
      }
    if ((ok == false) && (strcmp(argv[1], "-zeroRanksContribute") == 0))
      {
      argc--; argv++;
      ok = true;
      zeroRanksContribute = true;
      }
    if ((ok == false) && (strcmp(argv[1], "-r") == 0))
      {
      argc--; argv++;
      ok = true;
      fmRadius=atof(argv[1]);
      argc--; argv++;
      }

    if (ok == false) 
      {

      if (fileImage2D == NULL) 
        {
        fileImage2D = argv[1];
        argc--;
        argv++;
        }
      
      else if (fileVolume3D == NULL) 
        {
        fileVolume3D = argv[1];
        argc--;
        argv++;
        }
      
      else 
        {
        std::cerr << "ERROR: Can not parse argument " << argv[1] << std::endl;
        usage();
        }
      }
    } 

  if (verbose) 
    {
    if (fileImage2D)  std::cout << "Input 2D image: " << fileImage2D  << std::endl;
    if (fileVolume3D) std::cout << "Input 3D image: " << fileVolume3D << std::endl;
    if (fileOutput)   std::cout << "Output image: "   << fileOutput   << std::endl;
    if (fileOutputMasked)   std::cout << "Output image masked: "   << fileOutputMasked   << std::endl;
    }


  
  const    unsigned int    Dimension = 3;
  typedef  short           PixelType;
  typedef itk::Image< PixelType, Dimension > ImageType2D;
  typedef itk::Image< PixelType, Dimension > ImageType3D;
  typedef   float     InternalPixelType;
  typedef itk::Image< InternalPixelType, Dimension > InternalImageType;
  typedef   unsigned int     RankPixelType;
  typedef itk::Image< RankPixelType, Dimension > RankImageType;
  typedef itk::Euler3DTransform< double >  TransformType;
  typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
  typedef ora::StochasticRankCorrelationImageToImageMetric<InternalImageType,
        InternalImageType, RankImageType> MetricType;
  typedef itk::RayCastInterpolateImageFunction< 
                                    InternalImageType,
                                    double          >    InterpolatorType;
  typedef itk::MultiResolutionImageRegistrationMethod< 
                                    InternalImageType, 
                                    InternalImageType >   RegistrationType;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType > ImagePyramidType2D;
  typedef itk::MultiResolutionPyramidImageFilter<
                                    InternalImageType,
                                    InternalImageType > ImagePyramidType3D;

  ImagePyramidType2D::Pointer imagePyramid2D = 
      ImagePyramidType2D::New();
  ImagePyramidType3D::Pointer imagePyramid3D =
      ImagePyramidType3D::New();

  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();

  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  );

  registration->SetFixedImagePyramid( imagePyramid2D );
  registration->SetMovingImagePyramid( imagePyramid3D );
  
  if (debug) 
    {
    //metric->DebugOn();
    //transform->DebugOn();   
    optimizer->DebugOn();   
    //interpolator->DebugOn();
    //registration->DebugOn();
    }


  typedef itk::ImageFileReader< ImageType2D > ImageReaderType2D;
  typedef itk::ImageFileReader< ImageType3D > ImageReaderType3D;

  ImageReaderType2D::Pointer imageReader2D = ImageReaderType2D::New();
  ImageReaderType3D::Pointer imageReader3D = ImageReaderType3D::New();

  imageReader2D->SetFileName( fileImage2D );
  imageReader3D->SetFileName( fileVolume3D );

  typedef itk::CastImageFilter< 
                        ImageType2D, InternalImageType > CastFilterType2D;
  typedef itk::CastImageFilter< 
                        ImageType3D, InternalImageType > CastFilterType3D;

  CastFilterType2D::Pointer caster2D = CastFilterType2D::New();
  CastFilterType3D::Pointer caster3D = CastFilterType3D::New();
  caster2D->SetInput( imageReader2D->GetOutput() );
  caster3D->SetInput( imageReader3D->GetOutput() );
  caster2D->Update();
  caster3D->Update();
  registration->SetFixedImage(  caster2D->GetOutput() );
  registration->SetMovingImage( caster3D->GetOutput() );
  registration->SetFixedImageRegion( 
       caster2D->GetOutput()->GetBufferedRegion() );
   

  transform->SetComputeZYX(true);
  TransformType::OutputVectorType translation;
  translation[0] = tx;
  translation[1] = ty;
  translation[2] = tz;

  if (!initpars) // intuitive initialization
  {
    transform->SetTranslation(translation);
    transform->SetRotation(M_PI/180.0*rx, M_PI/180.0*ry, M_PI/180.0*rz);
  }
  else // set parameters directly
  {
    TransformType::ParametersType pp;
    pp.SetSize(6);
    for (int xx = 0; xx < 6; xx++)
      pp[xx] = initpars[xx];
    transform->SetParameters(pp);
  }

  double origin3D[ Dimension ];
  const itk::Vector<double, 3> resolution3D = imageReader3D->GetOutput()->GetSpacing();

  typedef ImageType3D::RegionType      ImageRegionType3D;
  typedef ImageRegionType3D::SizeType  SizeType3D;

  ImageRegionType3D region3D = caster3D->GetOutput()->GetBufferedRegion();
  SizeType3D        size3D   = region3D.GetSize();

  origin3D[0] = resolution3D[0]*((double) size3D[0])/2.; 
  origin3D[1] = resolution3D[1]*((double) size3D[1])/2.; 
  origin3D[2] = resolution3D[2]*((double) size3D[2])/2.; 

  TransformType::InputPointType center;
  center[0] = cx + origin3D[0];
  center[1] = cy + origin3D[1];
  center[2] = cz + origin3D[2];

  transform->SetCenter(center);
  // Software Guide : EndCodeSnippet

  if (verbose) 
    {
    std::cout << "3D image size: "
              << size3D[0] << ", " << size3D[1] << ", " << size3D[2] << std::endl
              << "   resolution: "
              << resolution3D[0] << ", " << resolution3D[1] << ", " << resolution3D[2] << std::endl
              << "Transform: " << transform << std::endl;
    }


  double origin2D[ Dimension ];
  const itk::Vector<double, 3> resolution2D = imageReader2D->GetOutput()->GetSpacing();

  typedef ImageType2D::RegionType      ImageRegionType2D;
  typedef ImageRegionType2D::SizeType  SizeType2D;

  ImageRegionType2D region2D = caster2D->GetOutput()->GetBufferedRegion();
  SizeType2D        size2D   = region2D.GetSize();

  origin2D[0] = origin3D[0] + o2Dx - resolution2D[0]*((double) size2D[0] - 1.)/2.; 
  origin2D[1] = origin3D[1] + o2Dy - resolution2D[1]*((double) size2D[1] - 1.)/2.; 
  origin2D[2] = origin3D[2] + focalLength/2.; 

  imageReader2D->GetOutput()->SetOrigin( origin2D );

  registration->SetFixedImageRegion( caster2D->GetOutput()->GetBufferedRegion() );

  if (verbose)
    {
    std::cout << "2D image size: "
              << size2D[0] << ", " << size2D[1] << ", " << size2D[2] << std::endl
              << "   resolution: "
              << resolution2D[0] << ", " << resolution2D[1] << ", " << resolution2D[2] << std::endl
              << "   and position: " 
              << origin2D[0] << ", " << origin2D[1] << ", " << origin2D[2] << std::endl;
    }


  std::vector<double> sampledRes2D;
  std::vector<ImageRegionType2D> imageRegionPyramid2D;

  imagePyramid2D->SetNumberOfLevels( nScales );
  
  sampledRes2D.resize( nScales ); // FIX by David Thivierge-gaulin - thank you!
  imageRegionPyramid2D.resize( nScales );

  typedef ImagePyramidType2D::ScheduleType   ScheduleType2D;
  ScheduleType2D schedule2D = imagePyramid2D->GetSchedule();

  for ( unsigned int dim = 0; dim < ImageType2D::ImageDimension; dim++ )
    {
    schedule2D[ nScales-1 ][ dim ] = maxScale;
    }

  for ( int level=nScales-2; level >= 0; level-- )
    {
    for ( unsigned int dim = 0; dim < ImageType2D::ImageDimension; dim++ )
      {
      schedule2D[ level ][ dim ] = 2*schedule2D[ level+1 ][ dim ];
      }
    }

  typedef ImageRegionType2D::IndexType       IndexType2D;
  IndexType2D       inputStart2D = region2D.GetIndex();

  for ( unsigned int level=0; level < nScales; level++ )
    {
    SizeType2D  size;
    IndexType2D start;
    sampledRes2D[ level ] = 0.;
    for ( unsigned int dim = 0; dim < ImageType2D::ImageDimension; dim++ )
      {
      float scaleFactor = static_cast<float>( schedule2D[ level ][ dim ] );

      size[ dim ] = static_cast<SizeType2D::SizeValueType>(
        floor( static_cast<float>( size2D[ dim ] ) / scaleFactor ) );

      if( size[ dim ] < 1 )
        {
        size[ dim ] = 1;
        schedule2D[ level ][ dim ] = 1;
        }

       std::cout << level << " " << dim << " " 
                << size[ dim ] << " " << schedule2D[ level ][ dim ] 
                << std::endl;

       scaleFactor = static_cast<float>( schedule2D[ level ][ dim ] );
       start[ dim ] = static_cast<IndexType2D::IndexValueType>(
         ceil(  static_cast<float>( inputStart2D[ dim ] ) / scaleFactor ) ); 

      if ( dim < 2 ) 
        {
        sampledRes2D[ level ] +=  scaleFactor*resolution2D[ dim ]
                                 *scaleFactor*resolution2D[ dim ];
        }
      }

    sampledRes2D[ level ] = sqrt( sampledRes2D[ level ] );

    imageRegionPyramid2D[ level ].SetSize( size );
    imageRegionPyramid2D[ level ].SetIndex( start );
    }

  imagePyramid2D->SetSchedule( schedule2D );


  std::vector<ImageRegionType3D> imageRegionPyramid3D;

  imagePyramid3D->SetNumberOfLevels( nScales );
  
  imageRegionPyramid3D.resize( nScales ); // FIX by David Thivierge-gaulin - thank you!

  typedef ImagePyramidType3D::ScheduleType   ScheduleType3D;
  ScheduleType3D schedule3D = imagePyramid3D->GetSchedule();
  
  for ( unsigned int level=0; level < nScales; level++ )
    {
    for ( unsigned int dim = 0; dim < ImageType3D::ImageDimension; dim++)
      {
      double res = resolution3D[ dim ];
      schedule3D[ level ][ dim ] = 1;
      while ( res*2. < sampledRes2D[ level ] ) 
        {
        schedule3D[ level ][ dim ] *= 2;
        res *= 2.;
        }
      }
    }
   

  typedef ImageRegionType3D::IndexType       IndexType3D;
  IndexType3D       inputStart3D = region3D.GetIndex();

  for ( unsigned int level=0; level < nScales; level++ )
    {
    SizeType3D  size;
    IndexType3D start;
    for ( unsigned int dim = 0; dim < ImageType3D::ImageDimension; dim++)
      {
      float scaleFactor = static_cast<float>( schedule3D[ level ][ dim ] );

      size[ dim ] = static_cast<SizeType3D::SizeValueType>(
        floor( static_cast<float>( size3D[ dim ] ) / scaleFactor ) );

      if( size[ dim ] < 1 )
        {
        size[ dim ] = 1;
        schedule3D[ level ][ dim ] = 1;
        }

      scaleFactor = static_cast<float>( schedule3D[ level ][ dim ] );
      start[ dim ] = static_cast<IndexType3D::IndexValueType>(
        ceil(  static_cast<float>( inputStart3D[ dim ] ) / scaleFactor ) ); 
      }
    imageRegionPyramid3D[ level ].SetSize( size );
    imageRegionPyramid3D[ level ].SetIndex( start );
    }

  imagePyramid3D->SetSchedule( schedule3D );


  InterpolatorType::InputPointType focalPoint;

  focalPoint[0] = origin3D[0];
  focalPoint[1] = origin3D[1];
  focalPoint[2] = origin3D[2] - focalLength/2.;

  interpolator->SetFocalPoint(focalPoint);
  interpolator->SetThreshold(threshold);
  interpolator->SetTransform(transform);

  if (verbose)
    {
    std::cout << "Focal point: " 
              << focalPoint[0] << ", " << focalPoint[1] << ", " << focalPoint[2] << std::endl
              << "Threshold: " << threshold << std::endl;
    }

  TransformType::ParametersType initialParameters = transform->GetParameters(); // store
  registration->SetInitialTransformParameters( transform->GetParameters() );
  optimizer->SetMaximize( false );
  optimizer->SetMaximumStepLength( maxStepSize );  
  optimizer->SetMinimumStepLength( minStepSize );
  optimizer->SetNumberOfIterations( maxit );
  optimizer->SetRelaxationFactor(rf);

  itk::Optimizer::ScalesType weightings( 6 );
  for (int xx = 0; xx < 6; xx++)
    weightings[xx] = oscale[xx];

  optimizer->SetScales( weightings );

  if (verbose)
    {
    optimizer->Print( std::cout );
    }

  MetricType::ScalesType mdscales;
  mdscales.SetSize(transform->GetNumberOfParameters());
  for (unsigned int d = 0; d < transform->GetNumberOfParameters(); d++)
    mdscales[d] = dscale[d];
  MetricType::SeedsType mrseeds;
  mrseeds.SetSize(MetricType::ImageDimension);
  for (unsigned int d = 0; d < MetricType::ImageDimension; d++)
    mrseeds[d] = rseeds[d];
  metric->SetRandomSeeds(mrseeds);
  metric->SetDerivativeScales(mdscales);
  metric->SetFixedHistogramClipAtEnds(fhistclip);
  metric->SetFixedHistogramMinIntensity(fhistmin);
  metric->SetFixedHistogramMaxIntensity(fhistmax);
  metric->SetFixedNumberOfHistogramBins(fhistbins);
  metric->SetMovingHistogramClipAtEnds(mhistclip);
  metric->SetMovingHistogramMinIntensity(mhistmin);
  metric->SetMovingHistogramMaxIntensity(mhistmax);
  metric->SetMovingNumberOfHistogramBins(mhistbins);
  metric->SetNoOverlapReactionMode(1);
  metric->SetNoOverlapMetricValue(1000);
  metric->SetSampleCoverage(samplecoverage);
  metric->SetUseHornTiedRanksCorrection(horn);
  metric->SetMovingZeroRanksContributeToMeasure(zeroRanksContribute);
  metric->SetExtractSampleDistribution(false);


  CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );

  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );
  registration->SetNumberOfLevels( nScales );

  imagePyramid3D->Print(std::cout);
  imagePyramid2D->Print(std::cout);

  // Set fixed image mask if radius is > 0
  typedef itk::Image<unsigned char,
      MetricType::FixedImageType::ImageDimension> MaskImageType;
  MaskImageType::Pointer fim = NULL;
  if (fmRadius > 0)
  {
    ImageType2D::ConstPointer fimage = imageReader2D->GetOutput();
    ImageType2D::RegionType freg =
        fimage->GetLargestPossibleRegion();

    fim = MaskImageType::New();
    fim->SetRegions(fimage->GetLargestPossibleRegion());
    fim->Allocate();
    fim->SetOrigin(fimage->GetOrigin());
    fim->SetSpacing(fimage->GetSpacing());
    fim->SetDirection(fimage->GetDirection());
    fim->FillBuffer(0);

    MaskImageType::SizeType fimsz = fim->GetLargestPossibleRegion().GetSize();
    typedef itk::ImageRegionIterator<MaskImageType> MaskIteratorType;
    MaskIteratorType fimit(fim, fim->GetLargestPossibleRegion());
    while (!fimit.IsAtEnd())
    {
      MaskImageType::IndexType idx = fimit.GetIndex();
      int xx = idx[0] - fimsz[0] / 2;
      int yy = idx[1] - fimsz[1] / 2;
      if ((xx * xx + yy * yy) <= (fmRadius * fmRadius))
        fimit.Set(1);
      else
        fimit.Set(0);
      ++fimit;
    }

    // set as fixed image mask
    typedef itk::ImageMaskSpatialObject<MaskImageType::ImageDimension>
        MaskSpatialObjectType;
    MaskSpatialObjectType::Pointer fspatial =
        MaskSpatialObjectType::New();
    fspatial->SetImage(fim);
    fspatial->Update();
    metric->SetFixedImageMask(fspatial);
  }


  if (verbose)
    {
    std::cout << "Starting the registration now..." << std::endl;
    }

  try 
    { 
    registration->StartRegistration(); 
    } 
  catch( itk::ExceptionObject & err ) 
    { 
      std::cout << "ExceptionObject caught !" << std::endl; 
      std::cout << err << std::endl; 
      return -1;
    } 

  typedef RegistrationType::ParametersType ParametersType;
  ParametersType finalParameters = registration->GetLastTransformParameters();

  const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
  const double bestValue = optimizer->GetValue();

  std::cout << "Result = " << std::endl;
  std::cout << " Rotation X = " << (finalParameters[0] / M_PI * 180) << " deg"  << std::endl;
  std::cout << " Rotation Y = " << (finalParameters[1] / M_PI * 180) << " deg"  << std::endl;
  std::cout << " Rotation Z = " << (finalParameters[2] / M_PI * 180) << " deg"  << std::endl;
  std::cout << " Translation X = " << finalParameters[3]  << std::endl;
  std::cout << " Translation Y = " << finalParameters[4]  << std::endl;
  std::cout << " Translation Z = " << finalParameters[5]  << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;

  
  typedef itk::ResampleImageFilter< InternalImageType, InternalImageType > FilterType;

  FilterType::Pointer filter = FilterType::New();

  filter->SetInput( caster3D->GetOutput() );
  filter->SetDefaultPixelValue( 0 );

  filter->SetTransform( transform );
  filter->SetInterpolator( interpolator );

  filter->SetSize( imageReader2D->GetOutput()->GetLargestPossibleRegion().GetSize() );
  filter->SetOutputOrigin(  imageReader2D->GetOutput()->GetOrigin() );
  filter->SetOutputSpacing( imageReader2D->GetOutput()->GetSpacing() );

  typedef itk::CastImageFilter< InternalImageType, ImageType2D > OutputCastFilterType;
  OutputCastFilterType::Pointer outputCaster   = OutputCastFilterType::New();
  outputCaster->SetInput( filter->GetOutput() );

  typedef itk::ImageFileWriter< ImageType2D >  WriterType;
  WriterType::Pointer writer = WriterType::New();

  writer->SetInput( outputCaster->GetOutput() );

  for (int xx = 0; xx < 3; xx++)
  {
    if (xx == 0 && !fileOutput)
      continue;

    if (xx == 1 && !initialOutput)
      continue;

    if (xx == 0)
    {
      writer->SetFileName( fileOutput );
    }
    else if (xx == 1)
    {
       writer->SetFileName( initialOutput );
       transform->SetParameters(initialParameters);
    }

    try { 
      std::cout << "Writing image: " << writer->GetFileName() << std::endl;
      writer->Update();
    } 
    catch( itk::ExceptionObject & err ) { 
      
      std::cerr << "ERROR: ExceptionObject caught !" << std::endl; 
      std::cerr << err << std::endl; 
    } 
  }

  if (fmRadius > 0)
  {
    // Write mask to file
    if (fileOutputMasked)
    {
      typedef itk::ImageFileWriter<MaskImageType> MaskWriterType;
      MaskWriterType::Pointer w = MaskWriterType::New();
      w->SetInput(fim);
      w->SetFileName("fixed_image_mask.mhd");
      try
      {
        w->Update();
      } catch (itk::ExceptionObject &e)
      {
        std::cerr << "Fixed image mask write error: " << e << "\n";
      }

       // Declare the type for the Mask image filter
       typedef itk::MaskImageFilter<ImageType2D, MaskImageType,
           ImageType2D > MaskFilterType;


       // Create a mask  Filter
       MaskFilterType::Pointer maskFilter = MaskFilterType::New();

       // Connect the input images
       maskFilter->SetInput1( outputCaster->GetOutput() );
       maskFilter->SetInput2( fim );
       maskFilter->SetOutsideValue( 0.0 );
       writer->SetInput(maskFilter->GetOutput());
       writer->SetFileName(fileOutputMasked);
       transform->SetParameters(finalParameters);
       try
       {
         writer->Update();
       } catch (itk::ExceptionObject &e)
       {
         std::cerr << "Fixed image mask write error: " << e << "\n";
       }
     }
    }

  if (initpars)
    delete[] initpars;

  return 0;
}

