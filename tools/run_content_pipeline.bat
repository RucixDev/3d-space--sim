@echo off
setlocal enabledelayedexpansion

set ROOT=%~dp0\..
set INPUT=%1
if "%INPUT%"=="" set INPUT=%ROOT%\tools\examples\universe_state.json
set OUTPUT=%2
if "%OUTPUT%"=="" set OUTPUT=%ROOT%\out\content
set SEED=%3
if "%SEED%"=="" set SEED=123

python -m tools.stellar_content_pipeline --input "%INPUT%" --output "%OUTPUT%" --seed %SEED%
echo Done. See: %OUTPUT%
