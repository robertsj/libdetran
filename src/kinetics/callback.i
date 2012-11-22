/* 
 *  Callback support largely unchanged from the contribution of "Carlos" at
 *    http://old.nabble.com/Automatically-generated-callbacks-using-directors-(python)-p19239301.html
 */

//----------------------------------------------------------------------------//
%define setCallbackMethod(num, klass, setter, oldArgs, newArgs, newParams, safe)
__setCallbackSupport__(num,                // unique id of macro instantiation
                       %extend klass {, }, // class containing callback setter
                       $self->setter,      // callback setter class method
                       setter,             // signature
                       oldArgs, 
                       newArgs, 
                       newParams, 
                       safe)
%enddef

//----------------------------------------------------------------------------//
%define __setCallbackSupport__(num, open, close, oldSetter, newSetter, oldArgs, 
                               newArgs, newParams, safe)

%feature("director") __Callback_##num##__;

%inline
%{
  class __Callback_##num##__ 
  {
  public:
    void safeCall##newArgs 
    {
      call##newParams;
    }
    virtual void call##newArgs = 0;
    virtual ~__Callback_##num##__() {} 
  };
%}

%{
  void __callbackWrapper_##num##__##oldArgs 
  {
    ((__Callback_##num##__*) data)->call##newParams;
  }
%}

%pythoncode 
%{
def __callable2Callback_##num##__(args):
  # len(args) == 2 -> we're in a method, else -> we're in a function
  if len(args) == 2: 
    arg = 1
  else: 
    arg = 0
  if not isinstance(args[arg], __Callback_##num##__):
    from operator import isCallable
    if not isCallable(args[arg]): 
      raise TypeError, "Neither callable nor Callback"
    callable = args[arg]
    class Callback_(__Callback_##num##__):
      def call(self, *args): callable(*args)
    if len(args) == 2: 
       args = (args[0], Callback_(),)
    else: 
      args = (Callback_(),)
  args[arg].__disown__()
  return args
%}

%pythonprepend newSetter %{args = __callable2Callback_##num##__(args)%}

open
    void newSetter(__Callback_##num##__* callback) 
    {
        oldSetter(__callbackWrapper_##num##__, (void*) callback);
    }
close

%enddef
