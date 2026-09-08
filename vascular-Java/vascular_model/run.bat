@echo off
chcp 65001 >nul
title Vascular Patterning Model

:: Check if compiled
if not exist "build\classes\vascular\VascularApp.class" (
    echo [ERROR] Not compiled yet. Run compile.bat first.
    pause
    exit /b 1
)

:menu
cls
echo ==============================================
echo   Vascular Patterning Model
echo   Turing-like Stochastic Reaction-Diffusion
echo   Based on Hearn (2019) PLoS ONE
echo ==============================================
echo.
echo   [1] GUI mode (interactive visualization)
echo   [2] Headless (default 2D HB model)
echo   [3] List all presets
echo   [4] 3D HBPM model (vascular bundles)
echo   [5] Large stem (diameter 80)
echo   [6] Headless + custom parameters
echo   [0] Exit
echo.
set /p choice="  Choice: "

if "%choice%"=="0" goto :end
if "%choice%"=="1" goto :gui
if "%choice%"=="2" goto :headless
if "%choice%"=="3" goto :list
if "%choice%"=="4" goto :hbpm
if "%choice%"=="5" goto :large
if "%choice%"=="6" goto :custom
goto :menu

:gui
echo Starting GUI mode...
java --enable-preview -cp build\classes vascular.VascularApp --gui
goto :menu

:headless
echo Running default 2D HB simulation...
java --enable-preview -cp build\classes vascular.VascularApp --headless
echo.
pause
goto :menu

:list
echo Available presets:
java --enable-preview -cp build\classes vascular.VascularApp --list
echo.
pause
goto :menu

:hbpm
echo Running 3D HBPM model (500M events, takes a while)...
java --enable-preview -cp build\classes vascular.VascularApp --preset 2
echo.
pause
goto :menu

:large
echo Running large stem (diameter 80)...
java --enable-preview -cp build\classes vascular.VascularApp --preset 4
echo.
pause
goto :menu

:custom
echo.
echo Enter parameter overrides (or press Enter to use defaults):
echo   Example: --DH=0.1 --k2=1.0 --DB=10.0
echo.
set /p params="  Parameters: "
echo.
echo Running with custom parameters: %params%
java --enable-preview -cp build\classes vascular.VascularApp --headless %params%
echo.
pause
goto :menu

:end
echo Goodbye!
