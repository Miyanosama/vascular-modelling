@rem Gradle wrapper script for Windows
@rem Generated manually for the Vascular Patterning Model project

@if "%DEBUG%"=="" @echo off
@rem Set local scope for the variables with windows NT shell
if "%OS%"=="Windows_NT" setlocal

set DIRNAME=%~dp0
if "%DIRNAME%"=="" set DIRNAME=.
@rem This is normally used to jump to the project root
set APP_BASE_NAME=%~n0
set APP_HOME=%DIRNAME%..

set WRAPPER_JAR=%APP_HOME%\gradle\wrapper\gradle-wrapper.jar

@rem Download gradle-wrapper.jar if not present
if not exist "%WRAPPER_JAR%" (
    echo.
    echo Downloading Gradle Wrapper JAR...
    echo.

    set DOWNLOAD_URL=https://raw.githubusercontent.com/gradle/gradle/v8.12.0/gradle/wrapper/gradle-wrapper.jar

    powershell -Command "& {
        $ProgressPreference = 'SilentlyContinue';
        [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12;
        Invoke-WebRequest -Uri '%DOWNLOAD_URL%' -OutFile '%WRAPPER_JAR%';
    }"

    if not exist "%WRAPPER_JAR%" (
        echo.
        echo ERROR: Could not download gradle-wrapper.jar automatically.
        echo Please download it manually from:
        echo   https://github.com/gradle/gradle/raw/v8.12.0/gradle/wrapper/gradle-wrapper.jar
        echo And place it at: %WRAPPER_JAR%
        echo.
        pause
        exit /b 1
    )
)

@rem Execute Gradle
set JAVA_EXE=javaw.exe
if defined JAVA_HOME set JAVA_EXE="%JAVA_HOME%\bin\javaw.exe"

%JAVA_EXE% ^
  "-Dorg.gradle.appname=%APP_BASE_NAME%" ^
  -classpath "%WRAPPER_JAR%" ^
  org.gradle.wrapper.GradleWrapperMain ^
  %*

@rem End local scope
if "%OS%"=="Windows_NT" endlocal

exit /b %ERRORLEVEL%
