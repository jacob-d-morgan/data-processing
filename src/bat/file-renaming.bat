@rem file-renaming.bat
@rem =================
@rem Recursively renames all the .caf and .did files in a directory.
@rem Prepends the date and time that they were created, and the name of the sub-directory they are in.


@echo off
setlocal enableDelayedExpansion

for %%F in ("C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Alan\*.did") do (
	set "MDate=%%~tF"
	set "ParsedDate=!MDate:~6,4!-!MDate:~3,2!-!MDate:~0,2!_!MDate:~11,2!!MDate:~14,2!"
	
	copy %%F C:\Users\Jacob\Documents\projects\data-processing\data\xp\copy-of-results\Results\Alan\Renamed\!ParsedDate!_%%~nxF )
rem exit