@rem prepend-date-created.bat
@rem ========================
@rem Recursively renames all the .caf and .did files in a folder.
@rem Prepends the date that they were created, and the name of the sub-directory they are in.
@rem ========================================================================================


rem Preamble
rem ========

@echo off
setlocal enableDelayedExpansion

cls
echo "File Dates" >fileDates.txt


rem Define Working Directory
rem ========================

set "rootDir=C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Test Folder\"

pause
cd %rootDir%

pause
mkdir "%rootDir%Processed-Files"

pause
set "getdate="


rem Loop Through Files
rem ==================

for /r "%rootDir%" %%F in ("*.did") do (
	
	rem Get Relative Path to File
	rem =========================
	
	set "absPath=%%~dpF
	for /l %%n in (1,1,500) do (if "!rootDir:~%%n,1!" neq "" (set /a "lenRootDir=%%n+1"))
	for /l %%N in (1,1,500) do (if "!absPath:~%%N,1!" neq "" (set /a "lenAbsPath=%%N"))
	for %%L in (!lenAbsPath!) do set "absPath=!absPath:~0,%%L!"
	for %%l in (!lenRootDir!) do set "relPath=!absPath:~%%l!"

	echo Root Dir: !rootDir! >>fileDates.txt
	echo Abs Path: !absPath! >>fileDates.txt
	echo Rel Path: !relPath! >>fileDates.txt

	rem Get Date File Created
	rem =====================

	set "filename=%%F"
	for /f "tokens=1-4delims=/:" %%a in ('dir /TC "!filename!"') do (
		if "%%c" neq "" (
			
			set "day=%%a"
			set "month=%%b"
			set "stringToParse=%%c"
			set "stringToParseAlso=%%d"
		
			set "year=!stringToParse:~0,4!"
			set "hour=!stringToParse:~6,2!"
			set "min=!stringToParseAlso:~0,2!"
			set "getdate=!year!-!month!-!day!_!hour!!min!"
			
			echo File Creation Date: !getdate! >>fileDates.txt
			echo.>>fileDates.txt
			echo.>>fileDates.txt))


	rem Generate New File Name and Copy to New Destination
	rem ==================================================
	
	set "folderNames=!relPath:\=--!
	set "folderNames=!folderNames: =-!

	set "newName=!getdate!_!folderNames!_%%~nxF
	
	
	copy "%%F" "!rootdir!Processed-Files\!newName!"

pause
	


echo Working in !relPath!...)

echo "Completed"