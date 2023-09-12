Echo OFF
REM dos script 20190813


set dmu1=dmu1.exe
set dmuai=dmuai.exe

set MKL_NUM_THREADS=2
IF "%1"=="-np" (
  set  =%2
  SHIFT & SHIFT
)

set FN=%1%.DIR
IF EXIST %FN% (  
  echo Starting DMU with %1.DIR as directive file
  IF %MKL_NUM_THREADS% GTR 1 (
    echo Running on %MKL_NUM_THREADS% CPU/Cores
  )
  %dmu1:"=% < %1.DIR > %1.lst
  IF EXIST MODINF (
    %dmuai:"=% >> %1.lst
  ) ELSE (
  echo 
  echo DMUAI Not started due to errors in DMU1 
   )
  IF EXIST %1.SOL  del %1.SOL
     rename SOL %1.SOL
  IF EXIST %1.PAROUT del %1.PAROUT
     rename PAROUT %1.PAROUT
  IF EXIST %1.PAROUT_STD del %1.PAROUT_std
     rename PAROUT_STD %1.PAROUT_STD
  IF EXIST %1.LLIK del %1.LLIK
     rename LLIK %1.LLIK
  IF EXIST RESIDUAL ( 
    IF EXIST %1.RESIDUAL del %1.RESIDUAL
       rename RESIDUAL %1.RESIDUAL
  )

  IF EXIST %1.INBREED del %1.INBREED

  IF EXIST INBREED ( 
    FOR %%I IN (INBREED) DO IF %%~zI equ 0 del INBREED ELSE rename INBREED %1.INBREED DONE
  )
  IF EXIST CODE_TABLE del CODE_TABLE
  IF EXIST DMU1.dir del DMU1.dir
  IF EXIST DMUAI.dir del DMUAI.dir
  IF EXIST DMU_LOG del DMU_LOG
  IF EXIST DUMMY del DUMMY
  IF EXIST FSPAKWK del FSPAKWK
  IF EXIST Latest_parm del Latest_parm
  IF EXIST LEVAL del LEVAL
  IF EXIST MODINF del MODINF
  IF EXIST PARIN del PARIN
  IF EXIST RCDATA_I del RCDATA_I
  IF EXIST RCDATA_R del RCDATA_R  
  FOR %%R IN (AINV*) DO del %%R DONE      
  FOR %%R IN (COR*) DO del %%R DONE
  FOR %%R IN (DOM*) DO del %%R DONE
  FOR %%R IN (IBD*) DO del %%R DONE
  FOR %%R IN (PEDFILE*) DO del %%R DONE
  FOR %%R IN (fort.*) DO del %%R DONE
 ) ELSE (
  ECHO ' '
  echo File %1.DIR not in current directory
)
