@echo off
setlocal enabledelayedexpansion

if "%1"=="" goto :usage
if not exist "%1" goto :no_file

set errors=0
set validated=0

for /f "tokens=1,* delims= " %%a in (%1) do call :check_file %%a "%%b"

echo.
echo Summary:
echo Files validated: %validated%
echo Files with errors: %errors%
goto :eof

:check_file
set hash=
for /f "skip=1 tokens=* usebackq" %%i in (`certutil -hashfile %2 MD5`) do (
  if not defined hash set "hash=%%i"
)
if not defined hash goto :hash_error
if "%hash%"=="%1" goto :file_ok
echo FAILED: %2
echo   Expected: %1
echo   Got: %hash%
set /a errors+=1
goto :eof

:hash_error
echo ERROR: Could not get hash for %2
set /a errors+=1
goto :eof

:file_ok
echo OK: %2
set /a validated+=1
goto :eof

:usage
echo Usage: validate_md5.bat md5summary.txt
goto :eof

:no_file
echo ERROR: MD5 summary file not found: %1
goto :eof
