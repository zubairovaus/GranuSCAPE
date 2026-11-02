@echo off
setlocal
set "GRANUSCAPE_PYTHON=%~dp0..\granuscape_env\Scripts\python.exe"
if not exist "%GRANUSCAPE_PYTHON%" (
  echo Run SETUP_TEST.ps1 first.
  pause
  exit /b 1
)
set /p "GRANUSCAPE_H5=Final HDF5 file or folder path, without quotes: "
set /p "GRANUSCAPE_QC_OUTPUT=New empty QC output folder path, without quotes: "
if not defined GRANUSCAPE_H5 exit /b 1
if not defined GRANUSCAPE_QC_OUTPUT exit /b 1
"%GRANUSCAPE_PYTHON%" -m granuscape qc --h5 "%GRANUSCAPE_H5%" --output "%GRANUSCAPE_QC_OUTPUT%" --config "%~dp0examples\qc.yaml"
set "GRANUSCAPE_EXIT=%ERRORLEVEL%"
pause
exit /b %GRANUSCAPE_EXIT%
