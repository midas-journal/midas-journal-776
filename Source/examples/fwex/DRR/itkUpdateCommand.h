
#ifndef __ITKUPDATECOMMAND_H
#define __ITKUPDATECOMMAND_H

#include "itkCommand.h"

namespace itk {

/**
 * \class UpdateCommand
 * \brief Updates a target at modified events.
 * \author Markus Neuner, neuner.markus e_T gmx.net
 *
 * This class receives Modified-events and calls Update on the specified
 * target at modified events (ModiefiedEvent).
 * The purpose of this class is to implement an automatic update of a class.
 *
 * \see Command
 * \see ModiefiedEvent
 */
template <class TTarget>
class UpdateCommand
  : public Command
{
public:
  /** Standard class typedefs. */
  typedef UpdateCommand       Self;
  typedef SmartPointer<Self>  Pointer;

  /** Typedefs of target. */
  typedef TTarget                      TargetType;
  typedef typename TargetType::Pointer TargetTypePointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(UpdateCommand, Command);

  /** Set target on what Update() gets called. */
  itkSetObjectMacro(Target, TargetType);

  /** Get target on what Update() gets called. */
  itkGetObjectMacro(Target, TargetType);

  /** Set whether automatic update is active or not. */
  itkSetMacro(UpdateActive, bool);

  /** Get whether automatic update is active or not. */
  itkGetConstMacro(UpdateActive, bool);

  /**  Update target if ModifiedEvent arrives and update is active. */
  virtual void Execute(Object *caller, const EventObject & event )
    {
    this->Execute((const itk::Object *) caller, event);
    }

  /**  Update target if ModifiedEvent arrives and update is active. */
  virtual void Execute( const Object *caller, const EventObject & event )
    {
    if (!m_UpdateActive || !ModifiedEvent().CheckEvent(&event))
      {
      return;
      }

    m_Target->Update();
    }

protected:
  /** Constructor: Initialise new instance. */
  UpdateCommand()
    {
    m_Target = NULL;
    m_UpdateActive = false;
    }

  virtual ~UpdateCommand(){}

  /** Class on what Update() gets called. **/
  TargetTypePointer m_Target;

  /** Determines whether automatic update is active or not. **/
  bool m_UpdateActive;

private:
  UpdateCommand(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk

#endif
