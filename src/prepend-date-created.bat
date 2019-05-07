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

set "dirToUse=C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Test Folder\Another Test Folder\Kenji 15 N"

set "getdate="


rem Loop Through Files
rem ==================

for %%F in ("%dirToUse%\*.did") do (
	
	rem Get Date File Created
	rem =====================

	set "filename=%%F"
	echo "File Name:">>fileDates.txt
	echo !filename!>>fileDates.txt
	
	echo "Output of 'dir'":>>fileDates.txt
	dir /TC "!filename!" >> fileDates.txt
	
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
			
			echo "File Creation Date:">>fileDates.txt
			echo !getdate!>>fileDates.txt
			echo.>>fileDates.txt))


	rem Get Relative Path to File
	rem =========================


echo "Working...")

echo "Completed"