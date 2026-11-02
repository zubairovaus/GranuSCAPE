@echo off
setlocal
set "GRANUSCAPE_PYTHON=%~dp0..\granuscape_env\Scripts\python.exe"
if not exist "%GRANUSCAPE_PYTHON%" (
  echo Run SETUP_TEST.ps1 first.
  pause
  exit /b 1
)
echo Uses the CZI folder selected in test.yaml. Full inference is repeated three times.
echo One additional first-image run is recorded separately as warm-up.
"%GRANUSCAPE_PYTHON%" -m granuscape benchmark --config "%~dp0test.yaml" --repeats 3 --warmup 1
set "GRANUSCAPE_EXIT=%ERRORLEVEL%"
pause
exit /b %GRANUSCAPE_EXIT%
