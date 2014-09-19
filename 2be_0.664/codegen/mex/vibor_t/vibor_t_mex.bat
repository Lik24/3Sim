@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2013a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2013a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=vibor_t_mex
set MEX_NAME=vibor_t_mex
set MEX_EXT=.mexw64
call mexopts.bat
echo # Make settings for vibor_t > vibor_t_mex.mki
echo COMPILER=%COMPILER%>> vibor_t_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> vibor_t_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> vibor_t_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> vibor_t_mex.mki
echo LINKER=%LINKER%>> vibor_t_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> vibor_t_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> vibor_t_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> vibor_t_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> vibor_t_mex.mki
echo BORLAND=%BORLAND%>> vibor_t_mex.mki
echo OMPFLAGS= >> vibor_t_mex.mki
echo OMPLINKFLAGS= >> vibor_t_mex.mki
echo EMC_COMPILER=msvcsdk>> vibor_t_mex.mki
echo EMC_CONFIG=optim>> vibor_t_mex.mki
"C:\Program Files\MATLAB\R2013a\bin\win64\gmake" -B -f vibor_t_mex.mk
