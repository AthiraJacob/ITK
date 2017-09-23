#include<itkImage.h>
#include<itkImageFileReader.h>
#include<itkImageFileWriter.h>
#include<itkMultiplyImageFilter.h>
#include<itkImageRegistrationMethod.h>
#include<itkMeanSquaresImageToImageMetric.h>
#include<itkLinearInterpolateImageFunction.h>
#include<itkRegularStepGradientDescentOptimizer.h>
#include<itkResampleImageFilter.h>
#include<itkAffineTransform.h>
#include<itkSubtractImageFilter.h>
#include "itkCastImageFilter.h"
#include<itkRescaleIntensityImageFilter.h>
#include<itkCommand.h>
#include<itkNearestNeighborInterpolateImageFunction.h>
#include "itkBSplineInterpolateImageFunction.h"

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

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    OptimizerPointer optimizer = 
                         dynamic_cast< OptimizerPointer >( object );

    if( ! itk::IterationEvent().CheckEvent( &event ) )
      {
      return;
      }

    std::cout << optimizer->GetCurrentIteration() << " = ";
    std::cout << optimizer->GetValue() << " : ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;
  }
   
};


int main(int argc, char *argv[])
{

typedef unsigned char pixtype;
const unsigned int dimension = 3;
typedef itk::Image<pixtype,dimension> InputImageType;
typedef itk::Image<pixtype,dimension> OutputImageType;

//read tr image & mask
typedef itk::ImageFileReader<InputImageType> ReaderType;
ReaderType::Pointer reader1=ReaderType::New();
reader1->SetFileName(argv[1]);
ReaderType::Pointer reader2=ReaderType::New();
reader2->SetFileName(argv[2]);
//multiply filter
typedef itk::MultiplyImageFilter<InputImageType,InputImageType> MultiplyFilter;  //just multiply?4 readers? output image pointer?
MultiplyFilter::Pointer multiplier1=MultiplyFilter::New();
multiplier1->SetInput1(reader1->GetOutput());
multiplier1->SetInput2(reader2->GetOutput());

//read head image & mask
ReaderType::Pointer reader3=ReaderType::New();
reader3->SetFileName(argv[3]);
ReaderType::Pointer reader4=ReaderType::New();
reader4->SetFileName(argv[4]);
//multiply filter
MultiplyFilter::Pointer multiplier2=MultiplyFilter::New();
multiplier2->SetInput1(reader3->GetOutput());
multiplier2->SetInput2(reader4->GetOutput());

//rescale
typedef itk::RescaleIntensityImageFilter<InputImageType,InputImageType >RescalerType;
//tr image
RescalerType::Pointer intensityRescaler1 = RescalerType::New();
intensityRescaler1->SetInput( multiplier1->GetOutput() );
intensityRescaler1->SetOutputMinimum( 0 );
intensityRescaler1->SetOutputMaximum( 255 );
intensityRescaler1->Update();
//head image
RescalerType::Pointer intensityRescaler2 = RescalerType::New();
intensityRescaler2->SetInput( multiplier2->GetOutput() );
intensityRescaler2->SetOutputMinimum( 0 );
intensityRescaler2->SetOutputMaximum( 255 );
intensityRescaler2->Update();


//registration
typedef itk::AffineTransform< double, dimension > TransformType;
typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;
typedef itk::MeanSquaresImageToImageMetric<InputImageType,InputImageType >    MetricType;
typedef itk::LinearInterpolateImageFunction<InputImageType, double >  InterpolatorType;
 typedef itk::ImageRegistrationMethod<InputImageType,InputImageType >RegistrationType;
 
  // Create components
  MetricType::Pointer         metric        = MetricType::New();
  TransformType::Pointer      transform     = TransformType::New();
  OptimizerType::Pointer      optimizer     = OptimizerType::New();
  InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
  RegistrationType::Pointer   registration  = RegistrationType::New();
 
  // Each component is now connected to the instance of the registration method.
  registration->SetMetric(        metric        );
  registration->SetOptimizer(     optimizer     );
  registration->SetTransform(     transform     );
  registration->SetInterpolator(  interpolator  ); 
  registration->SetFixedImage(intensityRescaler1->GetOutput());
  registration->SetMovingImage(intensityRescaler2->GetOutput());
  intensityRescaler1->Update();
  registration->SetFixedImageRegion(intensityRescaler1->GetOutput()->GetLargestPossibleRegion() );
  
  typedef RegistrationType::ParametersType ParametersType;
  ParametersType initialParameters( transform->GetNumberOfParameters() );
 
  // rotation matrix
  initialParameters[0] = 1.0;  
  initialParameters[1] = 0.0;  
  initialParameters[2] = 0.0;  
  initialParameters[3] = 0.0;  
  initialParameters[4] = 1.0;  // R(0,0)
  initialParameters[5] = 0.0;  // R(0,1)
  initialParameters[6] = 0.0;  // R(1,0)
  initialParameters[7] = 0.0;   
  initialParameters[8] = 1.0;
 
   
  // translation vector
  
  initialParameters[9] = 0.0;
  initialParameters[10] = 0.0;
  initialParameters[11] = 0.0;
 
  registration->SetInitialTransformParameters( initialParameters );
  
  optimizer->SetMaximumStepLength( 0.1 ); 
  optimizer->SetMinimumStepLength( 0.01 );
  optimizer->SetNumberOfIterations( 200 );
  
//scale parameters in optimizer
  typedef OptimizerType::ScalesType OptimizerScalesType;
OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );
const double translationScale = 1.0 / 2000.0;
optimizerScales[0] = 1.0;
optimizerScales[1] = 1.0;
optimizerScales[2] = 1.0;
optimizerScales[3] = 1.0;
optimizerScales[4] = 1.0;
optimizerScales[5] = 1.0;
optimizerScales[6] = 1.0;
optimizerScales[7] = 1.0;
optimizerScales[8] = 1.0;
optimizerScales[9] = translationScale;
optimizerScales[10] = translationScale;
optimizerScales[11] = translationScale;
optimizer->SetScales( optimizerScales );


CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
optimizer->AddObserver( itk::IterationEvent(), observer );
  
//register
   try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cerr << "ExceptionObject caught !" << std::endl;
    std::cerr << err << std::endl;
    return EXIT_FAILURE;
  }
  
ParametersType finalParameters = registration->GetLastTransformParameters();
 
  std::cout << "Final Parameters: " << finalParameters << std::endl;
 
  unsigned int numberOfIterations = optimizer->GetCurrentIteration();
 
  double bestValue = optimizer->GetValue();
 
  // Print out results
  std::cout << std::endl;
  std::cout << "Result = " << std::endl;
  std::cout << " Iterations    = " << numberOfIterations << std::endl;
  std::cout << " Metric value  = " << bestValue          << std::endl;
  
  
  //resampling moving image
typedef itk::ResampleImageFilter<InputImageType, InputImageType > ResampleFilterType;
ResampleFilterType::Pointer resampler = ResampleFilterType::New();

typedef itk::BSplineInterpolateImageFunction<InputImageType, double >  BsplineInterpolatorType;
BsplineInterpolatorType::Pointer Binterpolator = BsplineInterpolatorType::New();
resampler->SetInterpolator( Binterpolator );

resampler->SetInput(multiplier2->GetOutput() );
resampler->SetTransform(registration->GetOutput()->Get() );

InputImageType::Pointer fixedImage = multiplier1->GetOutput();
resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
resampler->SetOutputOrigin(fixedImage->GetOrigin() );
resampler->SetOutputSpacing(fixedImage->GetSpacing() );
resampler->SetOutputDirection(fixedImage->GetDirection() );
resampler->SetDefaultPixelValue(100);  
  
typedef itk::ImageFileWriter<InputImageType> WriterType;
WriterType::Pointer writer=WriterType::New();

//fixed image
writer->SetFileName("testtr.img");
writer->SetInput(intensityRescaler1->GetOutput());
writer->Update();

//registered image
writer->SetFileName("testreg.img");
writer->SetInput(resampler->GetOutput());
writer->Update();

//trimg-regimg
typedef itk::SubtractImageFilter< InputImageType, InputImageType, InputImageType > DifferenceFilterType;//SubtractImageFilter is used
DifferenceFilterType::Pointer difference = DifferenceFilterType::New();
difference->SetInput1( intensityRescaler1->GetOutput() );
difference->SetInput2( resampler->GetOutput() );
resampler->SetDefaultPixelValue( 1 );

//writing difference matrix
writer->SetFileName("testdiffb.img");
writer->SetInput(difference->GetOutput());
writer->Update();

RescalerType::Pointer intensityRescaler3 = RescalerType::New();
intensityRescaler3->SetInput( difference->GetOutput() );
intensityRescaler3->SetOutputMinimum(0 );
intensityRescaler3->SetOutputMaximum(255 );

//writing difference matrix
writer->SetFileName("testdiff.img");
writer->SetInput(intensityRescaler3->GetOutput());
writer->Update();

//transforming atlas
ReaderType::Pointer reader5=ReaderType::New();
reader5->SetFileName(argv[5]);

ResampleFilterType::Pointer atlasResampler = ResampleFilterType::New();

typedef itk::NearestNeighborInterpolateImageFunction<InputImageType, double >  atlasInterpolatorType;
atlasInterpolatorType::Pointer atlasInterpolator = atlasInterpolatorType::New();
atlasResampler->SetInterpolator( atlasInterpolator );
atlasResampler->SetInput(reader5->GetOutput() );
atlasResampler->SetTransform(registration->GetOutput()->Get() );
InputImageType::Pointer fixedImage1 = multiplier1->GetOutput();
atlasResampler->SetSize( fixedImage1->GetLargestPossibleRegion().GetSize() );
atlasResampler->SetOutputOrigin(fixedImage1->GetOrigin() );
atlasResampler->SetOutputSpacing(fixedImage1->GetSpacing() );
atlasResampler->SetOutputDirection(fixedImage1->GetDirection() );
atlasResampler->SetDefaultPixelValue(100);  

//writing transformed atlas
writer->SetFileName("testtrAtlas.img");
writer->SetInput(atlasResampler->GetOutput());
writer->Update();


return 0;
}

 
