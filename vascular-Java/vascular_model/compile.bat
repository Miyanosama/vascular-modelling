@echo off
chcp 65001 >nul
title Compile Vascular Model

echo ============================================
echo   Compiling Vascular Patterning Model
echo   JDK 25 + enable-preview
echo ============================================
echo.

set SRC=src\main\java\vascular
set OUT=build\classes

if not exist "%OUT%" mkdir "%OUT%"

echo [1/9] Compiling model package...
javac --enable-preview -source 25 -d %OUT% %SRC%\model\*.java
if errorlevel 1 goto :error

echo [2/9] Compiling reaction package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\reaction\*.java
if errorlevel 1 goto :error

echo [3/9] Compiling diffusion package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\diffusion\*.java
if errorlevel 1 goto :error

echo [4/9] Compiling simulation package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\simulation\*.java
if errorlevel 1 goto :error

echo [5/9] Compiling config package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\config\*.java
if errorlevel 1 goto :error

echo [6/9] Compiling visualization package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\visualization\*.java
if errorlevel 1 goto :error

echo [7/9] Compiling analysis package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\analysis\*.java
if errorlevel 1 goto :error

echo [8/9] Compiling io package...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\io\*.java
if errorlevel 1 goto :error

echo [9/9] Compiling main app...
javac --enable-preview -source 25 -d %OUT% -cp %OUT% %SRC%\VascularApp.java
if errorlevel 1 goto :error

echo.
echo ============================================
echo   COMPILE SUCCESSFUL!
echo   Classes written to: %OUT%
echo ============================================
echo.
echo Next: run run.bat to start the simulation
pause
goto :end

:error
echo.
echo ============================================
echo   COMPILE FAILED! See errors above.
echo ============================================
pause
exit /b 1

:end
