@echo off
setlocal
set "GRANUSCAPE_PYTHON=%~dp0..\granuscape_env\Scripts\python.exe"
if not exist "%GRANUSCAPE_PYTHON%" (
  echo Run SETUP_TEST.ps1 first.
  pause
  exit /b 1
)
"%GRANUSCAPE_PYTHON%" -m granuscape run --config "%~dp0test.yaml"
set "GRANUSCAPE_EXIT=%ERRORLEVEL%"
pause
exit /b %GRANUSCAPE_EXIT%
