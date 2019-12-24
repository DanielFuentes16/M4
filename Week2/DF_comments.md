# Comments lab 2

## Sift compilation and using
In order to use the sift library is necesary to change the paths where are the precompiled, then is necesary to compile again with mex.

- 1.  otool -L "".mexmaci64
  As you can see, it references MATLAB's libraries using @loader_path, which is wrong. That should be @rpath.
- 2. mex -compatibleArrayDims "".c
- 3.otool -L lbfgsC.mexmaci64
This looks a lot better, it's using @rpath as it should. The MEX-file now ran, meaning that the linker problem is solved.
  
  In some cases the .c is necessary to modify the import see.
  #include <bla.h>
is meant for standard library or framework headers, and the search strategy Is different than that used for

#include "bla.h"
See for example

What is the difference between #include <filename> and #include "filename"?
As a workaround, you can set the Xcode build setting "Always Search User Paths" to YES.