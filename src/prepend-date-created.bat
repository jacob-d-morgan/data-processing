@rem prepend-date-created.bat
@rem ========================
@rem Recursively renames all the .caf and .did files in a folder.
@rem Prepends the date that they were created, and the name of the sub-directory they are in.

cls

@echo off
setlocal enableDelayedExpansion
rem dir /tc "%filename%"

set "filename=C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Alan\ADRS_UZB2_57A_rep1-0115.did"

set "getdate="
for /f "tokens=1-4delims=/:" %%a in ('dir /TC C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Alan\*.did') do (
	if "%%c" neq "" (
		
		rem echo %%a
		rem echo %%b
		rem echo %%c
		rem echo %%d
		rem echo.

		set "day=%%a"
		set "month=%%b"
		set "stringToParse=%%c"
		set "stringToParseAlso=%%d"

		rem echo !stringToParse!
		
		set "year=!stringToParse:~0,4!"
		set "hour=!stringToParse:~6,2!"
		set "min=!stringToParseAlso:~0,2!"
		set "getdate=!year!-!month!-!day!_!hour!!min!" )

		echo.
		echo !getdate! )


rem echo "!getdate!"