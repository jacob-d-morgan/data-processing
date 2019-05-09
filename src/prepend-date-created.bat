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


rem Define Working Directory and Create Destination Directory
rem =========================================================

set "rootDir=C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Test Folder"
cd %rootDir%

for %%i in (.) do set "rootDirName=%%~nxi"
mkdir "%rootDir%\..\Processed-Files-from-%rootDirName%"

rem Loop Through Files
rem ==================

for /r "%rootDir%\" %%F in (*.caf *.did) do (

	rem Get Relative Path to File
	rem =========================
	
	rem Get the Path to the Directory from the File Path
	set "absPath=%%~dpF"

	rem Trim the trailing backslash from the Path to the Directory
	for /l %%N in (1,1,500) do (if "!absPath:~%%N,1!" neq "" (set /a "lenAbsPath=%%N"))
	for %%L in (!lenAbsPath!) do set "absPath=!absPath:~0,%%L!"
	
	rem Get the Relative Path from the Absolute Path
	for /l %%n in (1,1,500) do (if "!rootDir:~%%n,1!" neq "" (set /a "lenRootDir=%%n+2"))
	for %%l in (!lenRootDir!) do set "relPath=!absPath:~%%l!"


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
			)
		)


	rem Generate New File Name and Copy to New Destination
	rem ==================================================

	set "folderNames="
	if "!relPath!" neq "" (
		set "folderNames=!relPath:\=--!"
		set "folderNames=!folderNames: =-!"
		)

	set "newName=!getdate!_!folderNames!_%%~nxF"
	copy "%%F" "!rootDir!\..\Processed-Files-from-!rootDirName!\!newName!"
	
	echo Working: Copying files from !relPath!...
	)

echo "Script Completed"