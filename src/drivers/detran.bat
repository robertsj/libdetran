:: Utility script for setting Detran environment variables
:: on Windows systems.  Note, this does not permanantly change
:: the variables, but it does open a command window in which 
:: the variables are active.  Hence, from this window the 
:: user can execute python or xdetran.

set  PYTHONPATH=%PYTHONPATH%;%~dp0python
set  PATH=%PATH%;%~dp0bin;%~dp0lib;
call cmd