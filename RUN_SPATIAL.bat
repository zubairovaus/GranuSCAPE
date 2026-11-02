@echo off
setlocal
cd /d "%~dp0"
set "PYTHON=D:\sg_local\granuscape_env\Scripts\python.exe"
set "PYTHONPATH=%~dp0src;%PYTHONPATH%"
if not exist "%PYTHON%" (
 echo Python environment not found: %PYTHON%
 pause
 exit /b 1
)
set /p "SPATIAL_SOURCE=Paste final HDF5 folder path without quotes: "
set /p "SPATIAL_OUTPUT=Paste NEW output folder path without quotes: "
set /p "SPATIAL_SELECTION=Optional selection CSV path, or press Enter for all cells: "
if "%SPATIAL_SELECTION%"=="" (
 "%PYTHON%" -m granuscape spatial --h5 "%SPATIAL_SOURCE%" --output "%SPATIAL_OUTPUT%" --config "%~dp0examples\spatial.yaml"
) else (
 "%PYTHON%" -m granuscape spatial --h5 "%SPATIAL_SOURCE%" --output "%SPATIAL_OUTPUT%" --config "%~dp0examples\spatial.yaml" --selection "%SPATIAL_SELECTION%"
)
if errorlevel 1 echo ERROR: inspect the message above; outputs may be incomplete.
pause
