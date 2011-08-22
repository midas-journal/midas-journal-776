

#ifndef ITKEXHAUSTIVEOPTIMIZERCOMPLETE_HXX_
#define ITKEXHAUSTIVEOPTIMIZERCOMPLETE_HXX_


#include <itkExhaustiveOptimizer.h>


namespace itk
{


/** \class ExhaustiveOptimizerComplete
 * \brief This class extends the original ExhaustiveOptimizer by adding a getter
 * for the current iteration variable.
 * 
 * @see itk::ExhaustiveOptimizer
 * 
 * @author Philipp Steininger <phil.steininger e_T gmail.com>
 * @version 1.0
 */
class ITK_EXPORT ExhaustiveOptimizerComplete
  : public ExhaustiveOptimizer
{
public:
  /** Standard class typedefs. */
  typedef ExhaustiveOptimizerComplete Self;
  typedef ExhaustiveOptimizer Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;


  /** Run-time type information (and related methods). */
  itkTypeMacro(ExhaustiveOptimizerComplete, ExhaustiveOptimizer);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Getter for current iteration. **/
  itkGetConstReferenceMacro(CurrentIteration, unsigned long);

protected:
  /** Hidden constructor. **/
  ExhaustiveOptimizerComplete()
    : Superclass()
  {
  };

  /** Discarded destructor. **/
  virtual ~ExhaustiveOptimizerComplete() {};

  /** Print class information. **/
  virtual void PrintSelf(std::ostream& os, Indent indent) const
  {
    Superclass::PrintSelf(os, indent);

    os << indent << "ExhaustiveOptimizerComplete (extension)" << std::endl;
  };

private:
  ExhaustiveOptimizerComplete(const Self&); // purposely not implemented
  void operator=(const Self&); // purposely not implemented
};


}


#endif /* ITKEXHAUSTIVEOPTIMIZERCOMPLETE_HXX_ */
